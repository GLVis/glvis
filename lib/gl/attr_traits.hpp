#include "types.hpp"
#include "renderer_core.hpp"
#include <type_traits>

namespace gl3
{

struct attr_base
{
    static void setup() { }
    static void setup_legacy(void* buffer) { }
    static void clear() { }
    static void clear_legacy() { }
    enum { exists = false };
};

template<
    typename TV, typename TAttr, TAttr TV::*Attrib, typename TAttrInfo>
struct attr_exist
{
    constexpr static TAttr* get_attrib_offset()
    {
        return &(((TV*)0)->*Attrib);
    }
    // Sets up vertex attributes in the currently-bound buffer object.
    static void setup()
    {
        glEnableVertexAttribArray(TAttrInfo::AttrIdx);
        glVertexAttribPointer(TAttrInfo::AttrIdx,
                              std::tuple_size<TAttr>::value,
                              TAttrInfo::AttrGLType,
                              TAttrInfo::AttrNormalize,
                              sizeof(TV),
                              (void*) get_attrib_offset());
    }

    // Sets up client-side vertex pointers for the given buffer.
    static void setup_legacy(TV* buffer)
    {
        glEnableClientState(TAttrInfo::AttrGLArray);
        TAttrInfo::AttrFuncPtr(std::tuple_size<TAttr>::value,
                    TAttrInfo::AttrGLType,
                    sizeof(TV),
                    (char*) buffer + (size_t)get_attrib_offset());
    }

    static void clear()
    {
        glDisableVertexAttribArray(TAttrInfo::AttrIdx);
    }

    static void clear_legacy()
    {
        glDisableClientState(TAttrInfo::AttrGLArray);
    }
    
    enum { exists = true };
};

// Default attribute traits for vertex types
// Provides no-op setup/clear functions if an attribute doesn't exist.
template<typename TV, typename = int>
struct attr_coord : attr_base { };

template<typename TV, typename = int>
struct attr_normal : attr_base { };

template<typename TV, typename = int>
struct attr_color : attr_base { };

template<typename TV, typename = int>
struct attr_texcoord : attr_base { };

// Template specializations for attribute traits
// If an attribute exists in a vertex, generates setup and clear functions
// which setup OpenGL vertex attributes.
template<typename TV>
struct attr_coord<TV, decltype((void)TV::coord, 0)>
    : attr_exist<TV, decltype(TV::coord), &TV::coord,
                 attr_coord<TV, decltype((void)TV::coord, 0)>>
{
    const static int AttrIdx = CoreGLDevice::ATTR_VERTEX;
    const static bool AttrNormalize = false;

    const static GLenum AttrGLType = GL_FLOAT;
    const static GLenum AttrGLArray = GL_VERTEX_ARRAY;
    constexpr static auto AttrFuncPtr = &glVertexPointer;
};

template<typename TV>
struct attr_normal<TV, decltype((void)TV::norm, 0)>
    : attr_exist<TV, decltype(TV::norm), &TV::norm,
                 attr_normal<TV, decltype((void)TV::norm, 0)>>
{
    const static int AttrIdx = CoreGLDevice::ATTR_NORMAL;
    const static bool AttrNormalize = false;

    const static GLenum AttrGLType = GL_FLOAT;
    const static GLenum AttrGLArray = GL_NORMAL_ARRAY;
    static void AttrFuncPtr(GLint size, GLenum type, GLsizei stride, const GLvoid* ptr)
    {
        glNormalPointer(type, stride, ptr);
    }
};

template<typename TV>
struct attr_color<TV, decltype((void)TV::color, 0)>
    : attr_exist<TV, decltype(TV::color), &TV::color,
                 attr_color<TV, decltype((void)TV::color, 0)>>
{
    const static int AttrIdx = CoreGLDevice::ATTR_COLOR;
    const static bool AttrNormalize = true;

    const static GLenum AttrGLType = GL_UNSIGNED_BYTE;
    const static GLenum AttrGLArray = GL_COLOR_ARRAY;
    constexpr static auto AttrFuncPtr = &glColorPointer;
};

template<typename TV>
struct attr_texcoord<TV, decltype((void)TV::texCoord, 0)>
    : attr_exist<TV, decltype(TV::texCoord), &TV::texCoord,
                 attr_texcoord<TV, decltype((void)TV::texCoord, 0)>>
{
    const static int AttrIdx = CoreGLDevice::ATTR_TEXCOORD0;
    const static bool AttrNormalize = false;

    const static GLenum AttrGLType = GL_FLOAT;
    const static GLenum AttrGLArray = GL_TEXTURE_COORD_ARRAY;
    static void AttrFuncPtr(GLint size, GLenum type, GLsizei stride, const GLvoid* ptr)
    {
        glClientActiveTexture(GL_TEXTURE0);
        glTexCoordPointer(size, type, stride, ptr);
        glClientActiveTexture(GL_TEXTURE1);
        glTexCoordPointer(size, type, stride, ptr);
    }
};
}
