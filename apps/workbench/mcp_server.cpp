// ak_workbench — MCP server over HTTP on the loopback interface (see mcp_server.hpp).
#include "mcp_server.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <winsock2.h>
#include <ws2tcpip.h>
using socket_t = SOCKET;
using sock_len = int;
static void close_socket(socket_t s) { closesocket(s); }
static const socket_t kBad = INVALID_SOCKET;
#else
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <unistd.h>
using socket_t = int;
using sock_len = socklen_t;
static void close_socket(socket_t s) { ::close(s); }
static const socket_t kBad = -1;
#endif

namespace mcp {

namespace {

constexpr const char* kProtocol = "2025-06-18";

json rpc_result(const json& id, json result) { return {{"jsonrpc", "2.0"}, {"id", id}, {"result", std::move(result)}}; }

json rpc_error(const json& id, int code, const std::string& msg) {
    return {{"jsonrpc", "2.0"}, {"id", id}, {"error", {{"code", code}, {"message", msg}}}};
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

bool send_all(socket_t s, const std::string& data) {
    std::size_t off = 0;
    while (off < data.size()) {
        const int chunk = static_cast<int>(std::min<std::size_t>(data.size() - off, 1 << 20));
        const auto n = ::send(s, data.data() + off, chunk, 0);
        if (n <= 0) return false;
        off += static_cast<std::size_t>(n);
    }
    return true;
}

void reply(socket_t s, int status, const char* reason, const std::string& body,
           const char* type = "application/json") {
    std::string head = "HTTP/1.1 " + std::to_string(status) + " " + reason + "\r\n";
    if (!body.empty()) head += std::string("Content-Type: ") + type + "\r\n";
    head += "Content-Length: " + std::to_string(body.size()) + "\r\nConnection: close\r\n\r\n";
    send_all(s, head + body);
}

/// Origin allowed: absent (not a browser), or a localhost page.
bool origin_ok(const std::string& origin) {
    if (origin.empty() || origin == "null") return false;
    for (const char* ok : {"http://localhost", "https://localhost", "http://127.0.0.1", "https://127.0.0.1"})
        if (origin.rfind(ok, 0) == 0) {
            const char next = origin.size() > std::strlen(ok) ? origin[std::strlen(ok)] : '\0';
            if (next == '\0' || next == ':' || next == '/') return true;
        }
    return false;
}

} // namespace

json text(const std::string& s, bool error) {
    return {{"content", json::array({{{"type", "text"}, {"text", s}}})}, {"isError", error}};
}

json image_png(const std::vector<unsigned char>& png, const std::string& caption) {
    json content = json::array();
    content.push_back({{"type", "image"}, {"data", base64(png.data(), png.size())}, {"mimeType", "image/png"}});
    if (!caption.empty()) content.push_back({{"type", "text"}, {"text", caption}});
    return {{"content", content}, {"isError", false}};
}

std::string base64(const unsigned char* d, std::size_t n) {
    static const char* tbl = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve((n + 2) / 3 * 4);
    for (std::size_t i = 0; i < n; i += 3) {
        const unsigned v = (unsigned(d[i]) << 16) | (i + 1 < n ? unsigned(d[i + 1]) << 8 : 0u) |
                           (i + 2 < n ? unsigned(d[i + 2]) : 0u);
        out += tbl[(v >> 18) & 63];
        out += tbl[(v >> 12) & 63];
        out += i + 1 < n ? tbl[(v >> 6) & 63] : '=';
        out += i + 2 < n ? tbl[v & 63] : '=';
    }
    return out;
}

json Server::handle(const json& msg) {
    const bool is_request = msg.contains("id");
    const json id = is_request ? msg["id"] : json();
    if (!msg.is_object() || !msg.contains("method") || !msg["method"].is_string())
        return is_request ? rpc_error(id, -32600, "invalid request") : json();
    const std::string method = msg["method"];
    const json params = msg.value("params", json::object());
    if (!is_request) return json(); // notifications: initialized, cancelled, ...

    if (method == "initialize") {
        const std::string want = params.value("protocolVersion", std::string(kProtocol));
        return rpc_result(id, {{"protocolVersion", want.empty() ? kProtocol : want},
                               {"capabilities", {{"tools", {{"listChanged", false}}}}},
                               {"serverInfo", {{"name", "akernel-workbench"}, {"version", "0.1.0"}}},
                               {"instructions", instructions_}});
    }
    if (method == "ping") return rpc_result(id, json::object());
    if (method == "tools/list") {
        json list = json::array();
        for (const Tool& t : tools_)
            list.push_back({{"name", t.name}, {"description", t.description}, {"inputSchema", t.input_schema}});
        return rpc_result(id, {{"tools", list}});
    }
    if (method == "tools/call") {
        const std::string name = params.value("name", std::string());
        const json args = params.value("arguments", json::object());
        for (const Tool& t : tools_)
            if (t.name == name) {
                try {
                    return rpc_result(id, t.call(args));
                } catch (const std::exception& e) {
                    return rpc_result(id, text(std::string("error: ") + e.what(), true));
                }
            }
        return rpc_error(id, -32602, "unknown tool: " + name);
    }
    return rpc_error(id, -32601, "method not found: " + method);
}

bool Server::start(int port, std::string& error) {
    if (running_) return true;
#if defined(_WIN32)
    WSADATA wsa;
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) {
        error = "WSAStartup failed";
        return false;
    }
#endif
    socket_t s = ::socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (s == kBad) {
        error = "socket() failed";
        return false;
    }
    int yes = 1;
#if !defined(_WIN32)
    setsockopt(s, SOL_SOCKET, SO_REUSEADDR, reinterpret_cast<const char*>(&yes), static_cast<sock_len>(sizeof yes));
#else
    setsockopt(s, SOL_SOCKET, SO_EXCLUSIVEADDRUSE, reinterpret_cast<const char*>(&yes), static_cast<sock_len>(sizeof yes));
#endif
    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(static_cast<unsigned short>(port));
    addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    if (::bind(s, reinterpret_cast<sockaddr*>(&addr), static_cast<sock_len>(sizeof addr)) != 0 || ::listen(s, 8) != 0) {
        close_socket(s);
        error = "port " + std::to_string(port) + " is in use";
        return false;
    }
    listener_ = static_cast<std::intptr_t>(s);
    port_ = port;
    running_ = true;
    thread_ = std::thread([this] { loop(); });
    return true;
}

void Server::stop() {
    if (!running_) return;
    running_ = false;
    if (thread_.joinable()) thread_.join();
    close_socket(static_cast<socket_t>(listener_));
    listener_ = -1;
}

void Server::loop() {
    const socket_t ls = static_cast<socket_t>(listener_);
    while (running_) {
        fd_set set;
        FD_ZERO(&set);
        FD_SET(ls, &set);
        timeval tv{0, 200000};
        const int r = ::select(static_cast<int>(ls) + 1, &set, nullptr, nullptr, &tv);
        if (r <= 0) continue;
        sockaddr_in peer{};
        sock_len len = static_cast<sock_len>(sizeof peer);
        const socket_t c = ::accept(ls, reinterpret_cast<sockaddr*>(&peer), &len);
        if (c == kBad) continue;
        serve(static_cast<std::intptr_t>(c));
        close_socket(c);
    }
}

void Server::serve(std::intptr_t client) {
    const socket_t c = static_cast<socket_t>(client);
    // Read the header, then Content-Length bytes of body.
    std::string buf;
    char tmp[8192];
    std::size_t header_end = std::string::npos;
    while (header_end == std::string::npos) {
        const auto n = ::recv(c, tmp, static_cast<int>(sizeof tmp), 0);
        if (n <= 0) return;
        buf.append(tmp, static_cast<std::size_t>(n));
        header_end = buf.find("\r\n\r\n");
        if (buf.size() > (1u << 20) && header_end == std::string::npos) return;
    }
    const std::string head = buf.substr(0, header_end);
    std::string body = buf.substr(header_end + 4);
    const std::string request_line = head.substr(0, head.find("\r\n"));
    std::string method = request_line.substr(0, request_line.find(' '));
    std::string path = request_line.substr(method.size() + 1);
    path = path.substr(0, path.find(' '));
    std::size_t content_length = 0;
    std::string origin;
    std::size_t pos = head.find("\r\n");
    while (pos != std::string::npos && pos < head.size()) {
        const std::size_t next = head.find("\r\n", pos + 2);
        const std::string line = head.substr(pos + 2, (next == std::string::npos ? head.size() : next) - pos - 2);
        const std::size_t colon = line.find(':');
        if (colon != std::string::npos) {
            const std::string key = lower(line.substr(0, colon));
            std::string val = line.substr(colon + 1);
            val.erase(0, val.find_first_not_of(" \t"));
            if (key == "content-length") content_length = static_cast<std::size_t>(std::stoul(val));
            if (key == "origin") origin = val;
        }
        pos = next;
    }
    if (content_length > (64u << 20)) return reply(c, 413, "Payload Too Large", "");
    while (body.size() < content_length) {
        const auto n = ::recv(c, tmp, static_cast<int>(sizeof tmp), 0);
        if (n <= 0) return;
        body.append(tmp, static_cast<std::size_t>(n));
    }
    if (!origin.empty() && !origin_ok(origin)) return reply(c, 403, "Forbidden", "");
    if (path != "/mcp" && path != "/") return reply(c, 404, "Not Found", "");
    if (method == "GET") return reply(c, 405, "Method Not Allowed", ""); // no server-initiated stream
    if (method == "DELETE") return reply(c, 200, "OK", "");
    if (method != "POST") return reply(c, 405, "Method Not Allowed", "");

    ++requests_;
    json msg;
    try {
        msg = json::parse(body);
    } catch (const std::exception&) {
        return reply(c, 400, "Bad Request", rpc_error(nullptr, -32700, "parse error").dump());
    }
    if (msg.is_array()) { // a batch
        json out = json::array();
        for (const json& m : msg) {
            json r = handle(m);
            if (!r.is_null()) out.push_back(std::move(r));
        }
        if (out.empty()) return reply(c, 202, "Accepted", "");
        return reply(c, 200, "OK", out.dump());
    }
    json r = handle(msg);
    if (r.is_null()) return reply(c, 202, "Accepted", "");
    reply(c, 200, "OK", r.dump(-1, ' ', false, json::error_handler_t::replace));
}

} // namespace mcp
