#include <execinfo.h>
#include <fcntl.h>
#include <linux/fb.h>
#include <signal.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <unistd.h>

#include "glx.h"
#include <GLES/gl.h>

bool eglInitialized = false;
EGLDisplay eglDisplay;
EGLSurface eglSurface;
EGLConfig eglConfigs[1];

int8_t CheckEGLErrors() {
    EGLenum error;
    char *errortext;

    error = eglGetError();

    if (error != EGL_SUCCESS && error != 0)
    {
        switch (error)
        {
            case EGL_NOT_INITIALIZED:     errortext = "EGL_NOT_INITIALIZED"; break;
            case EGL_BAD_ACCESS:          errortext = "EGL_BAD_ACCESS"; break;
            case EGL_BAD_ALLOC:           errortext = "EGL_BAD_ALLOC"; break;
            case EGL_BAD_ATTRIBUTE:       errortext = "EGL_BAD_ATTRIBUTE"; break;
            case EGL_BAD_CONTEXT:         errortext = "EGL_BAD_CONTEXT"; break;
            case EGL_BAD_CONFIG:          errortext = "EGL_BAD_CONFIG"; break;
            case EGL_BAD_CURRENT_SURFACE: errortext = "EGL_BAD_CURRENT_SURFACE"; break;
            case EGL_BAD_DISPLAY:         errortext = "EGL_BAD_DISPLAY"; break;
            case EGL_BAD_SURFACE:         errortext = "EGL_BAD_SURFACE"; break;
            case EGL_BAD_MATCH:           errortext = "EGL_BAD_MATCH"; break;
            case EGL_BAD_PARAMETER:       errortext = "EGL_BAD_PARAMETER"; break;
            case EGL_BAD_NATIVE_PIXMAP:   errortext = "EGL_BAD_NATIVE_PIXMAP"; break;
            case EGL_BAD_NATIVE_WINDOW:   errortext = "EGL_BAD_NATIVE_WINDOW"; break;
            default:                      errortext = "unknown"; break;
        }

        printf("ERROR: EGL Error detected: %s (0x%X)\n", errortext, error);
        return 1;
    }

    return 0;
}

static int get_config_default(int attribute, int *value) {
    switch (attribute) {
        case GLX_USE_GL:
        case GLX_RGBA:
        case GLX_DOUBLEBUFFER:
            *value = 1;
            break;
        case GLX_STEREO:
            *value = 0;
            break;
        case GLX_AUX_BUFFERS:
            *value = 0;
            break;
        case GLX_RED_SIZE:
            *value = 5;
            break;
        case GLX_GREEN_SIZE:
            *value = 6;
            break;
        case GLX_BLUE_SIZE:
            *value = 5;
            break;
        case GLX_ALPHA_SIZE:
            *value = 8;
            break;
        case GLX_DEPTH_SIZE:
            *value = 16;
            break;
        case GLX_STENCIL_SIZE:
        case GLX_ACCUM_RED_SIZE:
        case GLX_ACCUM_GREEN_SIZE:
        case GLX_ACCUM_BLUE_SIZE:
        case GLX_ACCUM_ALPHA_SIZE:
            *value = 0;
            break;
        case GLX_RENDER_TYPE:
            *value = GLX_RGBA_BIT | GLX_COLOR_INDEX_BIT;
            break;
        case GLX_VISUAL_ID:
            *value = 1;
            break;
        case GLX_FBCONFIG_ID:
            *value = 1;
            break;
        case GLX_DRAWABLE_TYPE:
            *value = GLX_WINDOW_BIT;
            break;
        case 2: // apparently this is bpp
            *value = 16;
            return 0;
        default:
            printf("libGL: unknown attrib %i\n", attribute);
            *value = 0;
            return 1;
    }
    return 0;
}

// hmm...
EGLContext eglContext;
GLXContext glxContext;
Display *xDisplay;

#ifndef FBIO_WAITFORVSYNC
#define FBIO_WAITFORVSYNC _IOW('F', 0x20, __u32)
#endif
static bool g_showfps = false;
static bool g_usefb = false;
static bool g_vsync = false;
static bool g_xrefresh = false;
static bool g_stacktrace = false;
static int fbdev = -1;

static void xrefresh() {
    system("xrefresh");
}

static void signal_handler(int sig) {
    if (g_xrefresh)
        xrefresh();

    if (g_stacktrace) {
        switch (sig) {
            case SIGBUS:
            case SIGFPE:
            case SIGILL:
            case SIGSEGV: {
                void *array[10];
                size_t size = backtrace(array, 10);
                if (! size) {
                    printf("No stacktrace. Compile with -funwind-tables.\n");
                } else {
                    printf("Stacktrace: %i\n", size);
                    backtrace_symbols_fd(array, size, 2);
                }
                break;
            }
        }
    }
    signal(sig, SIG_DFL);
    raise(sig);
}

static void scan_env() {
    #define env(name, global, message)                    \
        char *env_##name = getenv(#name);                 \
        if (env_##name && strcmp(env_##name, "1") == 0) { \
            printf("libGL: " message "\n");               \
            global = true;                                \
        }

    env(LIBGL_XREFRESH, g_xrefresh, "xrefresh will be called on cleanup");
    env(LIBGL_STACKTRACE, g_stacktrace, "stacktrace will be printed on crash");
    if (g_xrefresh || g_stacktrace) {
        // TODO: a bit gross. Maybe look at this: http://stackoverflow.com/a/13290134/293352
        signal(SIGBUS, signal_handler);
        signal(SIGFPE, signal_handler);
        signal(SIGILL, signal_handler);
        signal(SIGSEGV, signal_handler);
        if (g_xrefresh) {
            signal(SIGINT, signal_handler);
            signal(SIGQUIT, signal_handler);
            signal(SIGTERM, signal_handler);
            atexit(xrefresh);
        }
    }
    env(LIBGL_FB, g_usefb, "framebuffer output enabled");
    env(LIBGL_FPS, g_showfps, "fps counter enabled");
    env(LIBGL_VSYNC, g_vsync, "vsync enabled");
    if (g_vsync) {
        fbdev = open("/dev/fb0", O_RDONLY);
        if (fbdev < 0) {
            fprintf(stderr, "Could not open /dev/fb0 for vsync.\n");
        }
    }
}

GLXContext glXCreateContext(Display *display,
                            XVisualInfo *visual,
                            GLXContext shareList,
                            Bool isDirect) {
    EGLint configAttribs[] = {
#ifdef PANDORA
        EGL_RED_SIZE, 5,
        EGL_GREEN_SIZE, 6,
        EGL_BLUE_SIZE, 5,
#endif
        EGL_DEPTH_SIZE, 16,
#ifdef USE_ES2
        EGL_RENDERABLE_TYPE, EGL_OPENGL_ES2_BIT,
#else
        EGL_BUFFER_SIZE, 16,
        EGL_SURFACE_TYPE, EGL_WINDOW_BIT,
        EGL_RENDERABLE_TYPE, EGL_OPENGL_ES_BIT,
#endif
        EGL_NONE
    };

#ifdef USE_ES2
    EGLint attrib_list[] = {
        EGL_CONTEXT_CLIENT_VERSION, 2,
        EGL_NONE
    };
#else
    EGLint *attrib_list = NULL;
#endif

    scan_env();
    GLXContext fake = malloc(sizeof(struct __GLXContextRec));
    if (eglDisplay != NULL) {
        eglMakeCurrent(eglDisplay, NULL, NULL, EGL_NO_CONTEXT);
        if (eglContext != NULL) {
            eglDestroyContext(eglDisplay, eglContext);
            eglContext = NULL;
        }
        if (eglSurface != NULL) {
            eglDestroySurface(eglDisplay, eglSurface);
            eglSurface = NULL;
        }
    }

    // make an egl context here...
    EGLBoolean result;
    if (eglDisplay == NULL || eglDisplay == EGL_NO_DISPLAY) {
        if (xDisplay == NULL) {
            xDisplay = display;
        }
        if (g_usefb) {
            eglDisplay = eglGetDisplay(EGL_DEFAULT_DISPLAY);
        } else {
            eglDisplay = eglGetDisplay(xDisplay);
        }
        if (eglDisplay == EGL_NO_DISPLAY) {
            printf("Unable to create EGL display.\n");
            return fake;
        }
    }

    // first time?
    if (eglInitialized == false) {
        eglBindAPI(EGL_OPENGL_ES_API);
        result = eglInitialize(eglDisplay, NULL, NULL);
        if (result != EGL_TRUE) {
            printf("Unable to initialize EGL display.\n");
            return fake;
        }
        eglInitialized = true;
    }

    int configsFound;
    result = eglChooseConfig(eglDisplay, configAttribs, eglConfigs, 1, &configsFound);
    CheckEGLErrors();
    if (result != EGL_TRUE || configsFound == 0) {
        printf("No EGL configs found.\n");
        return fake;
    }
    eglContext = eglCreateContext(eglDisplay, eglConfigs[0], EGL_NO_CONTEXT, attrib_list);
    CheckEGLErrors();

    // need to return a glx context pointing at it
    fake->display = xDisplay;
    fake->direct = true;
    fake->xid = 1;
    return fake;
}

void glXDestroyContext(Display *display, GLXContext ctx) {
    if (eglContext) {
        EGLBoolean result = eglDestroyContext(eglDisplay, eglContext);
        if (eglSurface != NULL) {
            eglDestroySurface(eglDisplay, eglSurface);
        }

        if (result != EGL_TRUE) {
            printf("Failed to destroy EGL context.\n");
        }
        if (fbdev >= 0) {
            close(fbdev);
            fbdev = -1;
        }
    }
    return;
}

Display *glXGetCurrentDisplay() {
    if (xDisplay && eglContext) {
        return xDisplay;
    }
    return NULL;
}

XVisualInfo *glXChooseVisual(Display *display,
                             int screen,
                             int *attributes) {

    // apparently can't trust the Display I'm passed?
    if (xDisplay == NULL) {
        xDisplay = display;
    }
    XVisualInfo *visual = (XVisualInfo *)malloc(sizeof(XVisualInfo));
    XMatchVisualInfo(display, screen, 16, TrueColor, visual);
    return visual;
}

/*
EGL_BAD_MATCH is generated if draw or read are not compatible with context
or if context is set to EGL_NO_CONTEXT and draw or read are not set to
EGL_NO_SURFACE, or if draw or read are set to EGL_NO_SURFACE and context is
not set to EGL_NO_CONTEXT.
*/

Bool glXMakeCurrent(Display *display,
                    int drawable,
                    GLXContext context) {

    if (eglDisplay != NULL) {
        eglMakeCurrent(eglDisplay, NULL, NULL, EGL_NO_CONTEXT);
        if (eglSurface != NULL) {
            eglDestroySurface(eglDisplay, eglSurface);
        }
    }

    if (g_usefb)
        drawable = 0;
    eglSurface = eglCreateWindowSurface(eglDisplay, eglConfigs[0], drawable, NULL);
    CheckEGLErrors();

    EGLBoolean result = eglMakeCurrent(eglDisplay, eglSurface, eglSurface, eglContext);
    CheckEGLErrors();
    if (result) {
        return true;
    }
    return false;
}

Bool glXMakeContextCurrent(Display *display, int drawable,
                           int readable, GLXContext context) {
    return glXMakeCurrent(display, drawable, context);
}

void glXSwapBuffers(Display *display,
                    int drawable) {
    render_raster();
    if (fbdev >= 0) {
        int arg = 0;
        ioctl(fbdev, FBIO_WAITFORVSYNC, &arg);
    }
    eglSwapBuffers(eglDisplay, eglSurface);
    CheckEGLErrors();

    if (g_showfps) {
        // framerate counter
        static int firstFrame = 0, lastFrame = 0;
        static int frames = 0;
        struct timeval out;
        gettimeofday(&out, NULL);
        frames++;
        if (!firstFrame) {
            lastFrame = firstFrame = out.tv_sec;
        } else if (out.tv_sec > lastFrame) {
            lastFrame = out.tv_sec;
            float fps = frames / (float)(lastFrame - firstFrame);
            printf("fps: %.2f\n", fps);
        }
    }
}

int glXGetConfig(Display *display,
                 XVisualInfo *visual,
                 int attribute,
                 int *value) {
    return get_config_default(attribute, value);
}

const char *glXQueryExtensionsString(Display *display, int screen) {
    const char *extensions = {
        "GLX_ARB_create_context "
        "GLX_ARB_create_context_profile "
        "GLX_EXT_create_context_es2_profile "
    };
    return extensions;
}

const char *glXQueryServerString(Display *display, int screen, int name) {
    return "";
}

Bool glXQueryExtension(Display *display, int *errorBase, int *eventBase) {
    if (errorBase)
        *errorBase = 0;

    if (eventBase)
        *eventBase = 0;

    return true;
}

Bool glXQueryVersion(Display *display, int *major, int *minor) {
    // TODO: figure out which version we want to pretend to implement
    *major = 1;
    *minor = 4;
    return true;
}

const char *glXGetClientString(Display *display, int name) {
    // TODO: return actual data here
    switch (name) {
        case GLX_VENDOR: break;
        case GLX_VERSION: break;
        case GLX_EXTENSIONS: break;
    }
    return "";
}

// stubs for glfw (GLX 1.3)
GLXContext glXGetCurrentContext() {
    // hack to make some games start
    return glxContext ? glxContext : (void *)1;
}

GLXFBConfig *glXChooseFBConfig(Display *display, int screen,
                       const int *attrib_list, int *count) {
    *count = 1;
    GLXFBConfig *configs = malloc(sizeof(GLXFBConfig) * *count);
    return configs;
}

GLXFBConfig *glXGetFBConfigs(Display *display, int screen, int *count) {
    *count = 1;
    GLXFBConfig *configs = malloc(sizeof(GLXFBConfig) * *count);
    return configs;
}

int glXGetFBConfigAttrib(Display *display, GLXFBConfig config, int attribute, int *value) {
    return get_config_default(attribute, value);
}

XVisualInfo *glXGetVisualFromFBConfig(Display *display, GLXFBConfig config) {
    if (xDisplay == NULL) {
        xDisplay = display;
    }
    XVisualInfo *visual = (XVisualInfo *)malloc(sizeof(XVisualInfo));
    XMatchVisualInfo(display, 0, 16, TrueColor, visual); 
    return visual;
}

GLXContext glXCreateNewContext(Display *display, GLXFBConfig config,
                               int render_type, GLXContext share_list,
                               Bool is_direct) {
    return glXCreateContext(display, 0, share_list, is_direct);
}

// misc stubs
void glXCopyContext(Display *display, GLXContext src, GLXContext dst, GLuint mask) {}
void glXCreateGLXPixmap(Display *display, XVisualInfo * visual, Pixmap pixmap) {} // should return GLXPixmap
void glXDestroyGLXPixmap(Display *display, void *pixmap) {} // really wants a GLXpixmap
void glXGetCurrentDrawable() {} // this should actually return GLXDrawable. Good luck.
Bool glXIsDirect(Display * display, GLXContext ctx) {
    return true;
}
void glXUseXFont(Font font, int first, int count, int listBase) {}
void glXWaitGL() {}
void glXWaitX() {}
