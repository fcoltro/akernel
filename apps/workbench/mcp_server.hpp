// ak_workbench — a small MCP server (Model Context Protocol, Streamable HTTP
// transport) so that an AI assistant such as Claude Code can drive the
// running workbench: create bodies, run the kernel's tools, read the audit
// results and take screenshots of the viewport.
//
// Transport: JSON-RPC 2.0 over HTTP POST to http://127.0.0.1:<port>/mcp,
// answered with a single application/json response (no SSE stream). The
// socket is bound to the loopback interface only, and requests carrying an
// Origin header from anywhere but localhost are refused (DNS rebinding).
//
// The server knows nothing about the kernel or the UI: main.cpp registers
// the tools. Tool handlers run on the server thread; handlers that touch the
// app's state hop to the UI thread themselves.
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#pragma warning(push, 0)
#endif
#include <json.hpp> // nlohmann/json, single header
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

namespace mcp {

using json = nlohmann::json;

struct Tool {
    std::string name;
    std::string description;
    json input_schema; ///< JSON Schema of the arguments (type: object)
    /// Returns the MCP CallToolResult: {"content": [...], "isError": bool}.
    std::function<json(const json& args)> call;
};

/// A text result, or an error result the assistant sees (tool errors are
/// results, not JSON-RPC errors, so the model can read and react to them).
json text(const std::string& s, bool error = false);
/// A PNG image result with a caption.
json image_png(const std::vector<unsigned char>& png, const std::string& caption);

class Server {
public:
    Server() = default;
    Server(const Server&) = delete;
    Server& operator=(const Server&) = delete;
    ~Server() { stop(); }

    void add(Tool t) { tools_.push_back(std::move(t)); }
    void set_instructions(std::string s) { instructions_ = std::move(s); }

    /// Listen on 127.0.0.1:port. False (with the reason) if the port is taken.
    bool start(int port, std::string& error);
    void stop();
    bool running() const { return running_; }
    int port() const { return port_; }
    int requests() const { return requests_; }

    /// Handle one JSON-RPC message; returns null for notifications.
    json handle(const json& msg);

private:
    void loop();
    void serve(std::intptr_t client);

    std::vector<Tool> tools_;
    std::string instructions_;
    std::thread thread_;
    std::atomic<bool> running_{false};
    std::atomic<int> requests_{0};
    std::intptr_t listener_ = -1;
    int port_ = 0;
};

/// Standard base64 (for image content).
std::string base64(const unsigned char* data, std::size_t n);

} // namespace mcp
