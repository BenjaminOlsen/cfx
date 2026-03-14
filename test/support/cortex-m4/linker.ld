/*
 * Linker script for Cortex-M4 bare-metal tests.
 * Target: LM3S6965EVB (qemu-system-arm -M lm3s6965evb)
 *
 * Memory layout:
 *   FLASH: 256K @ 0x00000000  (code + read-only data)
 *   RAM:    64K @ 0x20000000  (stack, .data, .bss, heap)
 */

MEMORY {
    FLASH (rx)  : ORIGIN = 0x00000000, LENGTH = 256K
    RAM   (rwx) : ORIGIN = 0x20000000, LENGTH = 64K
}

_estack = ORIGIN(RAM) + LENGTH(RAM);

SECTIONS {
    .isr_vector : {
        KEEP(*(.isr_vector))
    } > FLASH

    .text : {
        *(.text*)
    } > FLASH

    .rodata : {
        *(.rodata*)
    } > FLASH

    _sidata = LOADADDR(.data);

    .data : {
        _sdata = .;
        *(.data*)
        _edata = .;
    } > RAM AT > FLASH

    .bss : {
        _sbss = .;
        *(.bss*)
        *(COMMON)
        _ebss = .;
    } > RAM

    _end = .;
}
