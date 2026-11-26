// The stb_image_write implementation, in its own target so the kernel's
// strict warning flags do not apply to third-party code.
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_WRITE_NO_STDIO
#include <stb_image_write.h>
