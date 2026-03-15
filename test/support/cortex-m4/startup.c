/*
 * Cortex-M4 bare-metal startup code for QEMU semihosting tests.
 *
 * Provides vector table, .data/.bss init, and semihosting exit.
 * Target board: LM3S6965EVB (qemu-system-arm -M lm3s6965evb)
 */
#include <stdint.h>

extern uint32_t _estack;
extern uint32_t _sdata, _edata, _sidata;
extern uint32_t _sbss, _ebss;
extern int main(void);

void Reset_Handler(void)
{
    /* Copy .data from flash to RAM */
    uint32_t *src = &_sidata;
    uint32_t *dst = &_sdata;
    while (dst < &_edata) {
        *dst++ = *src++;
    }

    /* Zero .bss */
    dst = &_sbss;
    while (dst < &_ebss) {
        *dst++ = 0;
    }

    int ret = main();

    /*
     * Exit via ARM semihosting.
     * r0 = 0x18 (SYS_EXIT / angel_SWIreason_ReportException)
     * r1 = exit reason (0x20026 = ADP_Stopped_ApplicationExit on success,
     *                   0x20000 = ADP_Stopped_InternalError on failure)
     *
     * Note: can't use "i" constraint here because ret is a runtime value.
     * Instead, compute the semihosting exit code in a register.
     */
    uint32_t exit_code = (ret == 0) ? 0x20026u : 0x20000u;
    __asm__ volatile (
        "mov r0, #0x18\n"
        "mov r1, %0\n"
        "bkpt #0xAB\n"
        : : "r" (exit_code) : "r0", "r1"
    );
    while (1) {}
}

/* Minimal vector table - only reset handler needed for testing */
__attribute__((section(".isr_vector")))
void (* const vectors[])(void) = {
    (void (*)(void))&_estack,   /* Initial stack pointer */
    Reset_Handler,               /* Reset handler */
};
