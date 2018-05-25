#ifndef SDL_HPP
#define SDL_HPP
#include <string>
#include <memory>
#include <functional>
#include <map>
#include <SDL2/SDL.h>

struct EventInfo {
    GLint event;
    GLint data[4];
};

class SdlWindow {

public:
    typedef std::function<void(EventInfo*)> MouseDelegate;
    typedef std::function<void(GLenum)> KeyDelegate;
    typedef std::function<void(int, int)> WindowDelegate;
    typedef std::function<void()> Delegate;

    static KeyDelegate toKeyFn(Delegate func) {
        return [func](GLenum e) { func(); }
    }
private:
    struct _SdlHandle;
    //Use shared_ptr to manage handle lifetimes
    std::shared_ptr<_SdlHandle> _handle;
    Delegate onIdle;
    WindowDelegate onExpose
    WindowDelegate onReshape;
    std::map<int, KeyDelegate> onKeyDown;
    std::map<int, MouseDelegate> onMouseDown;
    std::map<int, MouseDelegate> onMouseUp;
    std::map<int, MouseDelegate> onMouseMove;

public:
    SdlWindow(std::string& title, int w, int h);
    ~SdlWindow();

    /**
     * Creates a new OpenGL context with the window.
     * Returns false if OpenGL or GLEW intialization fails.
     */
    bool createGlContext();
    /**
     * Runs the window loop.
     */
    void mainLoop();
    void setOnIdle(Delegate func) { onIdle = func; }
    void setOnExpose(WindowDelegate func) { onExpose = func; }
    void setOnReshape(WindowDelegate func) { onReshape = func; }
    
    void setOnKeyDown(SDL_Keycode key, Delegate func) {
        onKeyDown[key] = [func](GLenum e) { func(); }
    }
    void setOnKeyDown(SDL_Keycode key, KeyDelegate upper, KeyDelegate lower) {
        onKeyDown[key] = [upper, lower](GLenum e) {
            if ((e & KMOD_SHIFT) && upper) { upper(e); }
            else { if (lower) lower(e); }
        }
    }
    void setOnKeyDown(SDL_Keycode key, Delegate upper, Delegate lower) {
        onKeyDown[key] = [upper, lower](GLenum e) {
            if ((e & KMOD_SHIFT) && upper) { upper(); }
            else { if (lower) lower(); }
        }
    }
    void setOnKeyDown(SDL_Keycode key, KeyDelegate func) { onKeyDown[key] = func; }
    
    void setOnMouseDown(int btn, MouseDelegate func) { onMouseDown[btn] = func; }
    void setOnMouseUp(int btn, MouseDelegate func) { onMouseUp[btn] = func; }
    void setOnMouseMove(int btn, MouseDelegate func) { onMouseMove[btn] = func; }
    
    void callKeyDown(SDL_Keycode k) { onKeyDown[k](0); } 

    void setWindowTitle(std::string& title);
    void setWindowSize(int w, int h);
    void setWindowPos(int x, int y);

    operator bool() { return _handle; }
    bool isWindowInitialized() { return _handle; }
    /**
     * Returns true if the OpenGL context was successfully initialized.
     */
    bool isGlInitialized();
};

#endif
