/*
 * Semihosting syscall stubs for Cortex-M4 bare-metal tests.
 *
 * These replace the newlib nosys.specs stubs with semihosting calls,
 * enabling printf output and proper exit codes under QEMU.
 */
#include <stdint.h>
#include <errno.h>
#include <sys/stat.h>

/* Semihosting operation numbers */
#define SYS_OPEN   0x01
#define SYS_CLOSE  0x02
#define SYS_WRITEC 0x03
#define SYS_WRITE  0x05
#define SYS_READ   0x06
#define SYS_EXIT   0x18

static inline int semihost_call(int op, void *arg)
{
    register int r0 __asm__("r0") = op;
    register void *r1 __asm__("r1") = arg;
    __asm__ volatile ("bkpt #0xAB" : "+r"(r0) : "r"(r1) : "memory");
    return r0;
}

/*
 * Semihosting file handles for stdout/stderr.
 * Must be opened via SYS_OPEN(":tt", mode) — raw Unix fd numbers don't work.
 *   mode 4 = "w" (stdout), mode 8 = "a" (stderr)
 */
static int sh_stdout = -1;
static int sh_stderr = -1;

static int sh_open_tt(int mode)
{
    uint32_t args[3] = { (uint32_t)(uintptr_t)":tt", (uint32_t)mode, 3 /* strlen(":tt") */ };
    return semihost_call(SYS_OPEN, args);
}

static void sh_init_handles(void)
{
    if (sh_stdout < 0) sh_stdout = sh_open_tt(4);  /* "w" */
    if (sh_stderr < 0) sh_stderr = sh_open_tt(8);  /* "a" */
}

int _write(int fd, const void *buf, int len)
{
    sh_init_handles();
    int handle = (fd == 2) ? sh_stderr : sh_stdout;
    uint32_t args[3] = { (uint32_t)handle, (uint32_t)buf, (uint32_t)len };
    int remaining = semihost_call(SYS_WRITE, args);
    return len - remaining;
}

void _exit(int status)
{
    uint32_t code = (status == 0) ? 0x20026u : 0x20000u;
    semihost_call(SYS_EXIT, (void *)code);
    while (1) {}
}

void *_sbrk(int incr)
{
    extern char _end;
    static char *heap_end = 0;
    if (heap_end == 0)
        heap_end = &_end;

    /* Guard: don't let heap grow into the stack (256-byte safety margin) */
    register char *sp __asm__("sp");
    if (heap_end + incr > sp - 256) {
        errno = ENOMEM;
        return (void *)-1;
    }

    char *prev = heap_end;
    heap_end += incr;
    return prev;
}

int _close(int fd)                          { (void)fd; return -1; }
int _fstat(int fd, struct stat *st)         { (void)fd; st->st_mode = S_IFCHR; return 0; }
int _isatty(int fd)                         { (void)fd; return 1; }
int _lseek(int fd, int ptr, int dir)        { (void)fd; (void)ptr; (void)dir; return 0; }
int _read(int fd, void *buf, int len)       { (void)fd; (void)buf; (void)len; return 0; }
int _kill(int pid, int sig)                 { (void)pid; (void)sig; return -1; }
int _getpid(void)                           { return 1; }
int _open(const char *name, int flags, int mode) { (void)name; (void)flags; (void)mode; return -1; }

#include <sys/times.h>
clock_t _times(struct tms *buf)             { if (buf) { buf->tms_utime = 0; buf->tms_stime = 0; buf->tms_cutime = 0; buf->tms_cstime = 0; } return 0; }
