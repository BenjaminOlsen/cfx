	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 13, 0	sdk_version 13, 3
	.intel_syntax noprefix
	.globl	_cfx_poly1305_init              ## -- Begin function cfx_poly1305_init
	.p2align	4, 0x90
_cfx_poly1305_init:                     ## @cfx_poly1305_init
	.cfi_startproc
## %bb.0:
	push	rbp
	.cfi_def_cfa_offset 16
	.cfi_offset rbp, -16
	mov	rbp, rsp
	.cfi_def_cfa_register rbp
	mov	eax, 67108863
	and	eax, dword ptr [rsi]
	mov	dword ptr [rdi], eax
	mov	eax, dword ptr [rsi + 3]
	shr	eax, 2
	and	eax, 67108611
	mov	dword ptr [rdi + 4], eax
	mov	ecx, dword ptr [rsi + 6]
	shr	ecx, 4
	and	ecx, 67092735
	mov	dword ptr [rdi + 8], ecx
	mov	edx, dword ptr [rsi + 9]
	shr	edx, 6
	and	edx, 66076671
	mov	dword ptr [rdi + 12], edx
	mov	r8d, dword ptr [rsi + 12]
	shr	r8d, 8
	and	r8d, 1048575
	mov	dword ptr [rdi + 16], r8d
	lea	rax, [rax + 4*rax]
	mov	qword ptr [rdi + 24], rax
	lea	rax, [rcx + 4*rcx]
	mov	qword ptr [rdi + 32], rax
	lea	rax, [rdx + 4*rdx]
	mov	qword ptr [rdi + 40], rax
	lea	rax, [r8 + 4*r8]
	mov	qword ptr [rdi + 48], rax
	mov	eax, dword ptr [rsi + 16]
	mov	dword ptr [rdi + 56], eax
	mov	eax, dword ptr [rsi + 20]
	mov	dword ptr [rdi + 60], eax
	mov	eax, dword ptr [rsi + 24]
	mov	dword ptr [rdi + 64], eax
	mov	eax, dword ptr [rsi + 28]
	mov	dword ptr [rdi + 68], eax
	vxorps	xmm0, xmm0, xmm0
	vmovups	xmmword ptr [rdi + 72], xmm0
	mov	dword ptr [rdi + 88], 0
	pop	rbp
	ret
	.cfi_endproc
                                        ## -- End function
	.globl	_cfx_poly1305_update            ## -- Begin function cfx_poly1305_update
	.p2align	4, 0x90
_cfx_poly1305_update:                   ## @cfx_poly1305_update
	.cfi_startproc
## %bb.0:
	push	rbp
	.cfi_def_cfa_offset 16
	.cfi_offset rbp, -16
	mov	rbp, rsp
	.cfi_def_cfa_register rbp
	push	r15
	push	r14
	push	r13
	push	r12
	push	rbx
	sub	rsp, 120
	.cfi_offset rbx, -56
	.cfi_offset r12, -48
	.cfi_offset r13, -40
	.cfi_offset r14, -32
	.cfi_offset r15, -24
	mov	rcx, rdx
	mov	rax, qword ptr [rip + ___stack_chk_guard@GOTPCREL]
	mov	rax, qword ptr [rax]
	mov	qword ptr [rbp - 48], rax
	mov	r14d, dword ptr [rdi]
	mov	eax, dword ptr [rdi + 4]
	mov	qword ptr [rbp - 104], rax      ## 8-byte Spill
	mov	eax, dword ptr [rdi + 8]
	mov	qword ptr [rbp - 128], rax      ## 8-byte Spill
	mov	eax, dword ptr [rdi + 12]
	mov	qword ptr [rbp - 120], rax      ## 8-byte Spill
	mov	eax, dword ptr [rdi + 16]
	mov	qword ptr [rbp - 112], rax      ## 8-byte Spill
	mov	rdx, qword ptr [rdi + 24]
	mov	rax, qword ptr [rdi + 32]
	mov	qword ptr [rbp - 80], rax       ## 8-byte Spill
	mov	r8, qword ptr [rdi + 40]
	mov	r10, qword ptr [rdi + 48]
	mov	r11d, dword ptr [rdi + 72]
	mov	eax, dword ptr [rdi + 76]
	mov	r12d, dword ptr [rdi + 80]
	mov	r13d, dword ptr [rdi + 84]
	mov	qword ptr [rbp - 136], rdi      ## 8-byte Spill
	mov	r15d, dword ptr [rdi + 88]
	mov	rdi, rcx
	cmp	rcx, 16
	mov	qword ptr [rbp - 88], rdx       ## 8-byte Spill
	mov	qword ptr [rbp - 72], r14       ## 8-byte Spill
	jb	LBB1_2
	.p2align	4, 0x90
LBB1_1:                                 ## =>This Inner Loop Header: Depth=1
	mov	qword ptr [rbp - 152], rsi      ## 8-byte Spill
	mov	qword ptr [rbp - 96], rdi       ## 8-byte Spill
	mov	r9d, dword ptr [rsi]
	mov	ecx, 67108863
	and	r9d, ecx
	add	r9d, r11d
	mov	edi, dword ptr [rsi + 3]
	mov	ebx, dword ptr [rsi + 6]
	shr	edi, 2
	and	edi, 67108863
	add	edi, eax
	shr	ebx, 4
	and	ebx, 67108863
	add	ebx, r12d
	mov	r11d, dword ptr [rsi + 9]
	shr	r11d, 6
	add	r11d, r13d
	mov	eax, dword ptr [rsi + 12]
	shr	eax, 8
	mov	rsi, rdx
	lea	edx, [r15 + rax]
	add	edx, 16777216
	mov	rcx, r9
	imul	rcx, r14
	mov	r14, r10
	imul	r14, rdi
	add	r14, rcx
	mov	rcx, r8
	imul	rcx, rbx
	mov	r12, qword ptr [rbp - 80]       ## 8-byte Reload
	mov	r15, r12
	imul	r15, r11
	add	r15, rcx
	add	r15, r14
	imul	rsi, rdx
	add	rsi, r15
	mov	qword ptr [rbp - 144], rsi      ## 8-byte Spill
	mov	rcx, r9
	mov	rax, qword ptr [rbp - 104]      ## 8-byte Reload
	imul	rcx, rax
	mov	r15, rdi
	imul	r15, qword ptr [rbp - 72]       ## 8-byte Folded Reload
	add	r15, rcx
	mov	rcx, r10
	imul	rcx, rbx
	mov	r13, r8
	imul	r13, r11
	add	r13, rcx
	add	r13, r15
	mov	r15, r12
	imul	r15, rdx
	add	r15, r13
	mov	rcx, r9
	mov	rsi, qword ptr [rbp - 128]      ## 8-byte Reload
	imul	rcx, rsi
	mov	r13, rdi
	mov	r12, rax
	imul	r13, rax
	add	r13, rcx
	mov	rcx, rbx
	imul	rcx, qword ptr [rbp - 72]       ## 8-byte Folded Reload
	mov	r14, r10
	mov	r10, r8
	mov	r8, r14
	imul	r8, r11
	add	r8, rcx
	add	r8, r13
	mov	r13, r10
	imul	r13, rdx
	add	r13, r8
	mov	rcx, r9
	mov	rax, qword ptr [rbp - 120]      ## 8-byte Reload
	imul	rcx, rax
	mov	r8, rdi
	imul	r8, rsi
	add	r8, rcx
	mov	rcx, rbx
	imul	rcx, r12
	mov	r12, r11
	imul	r12, qword ptr [rbp - 72]       ## 8-byte Folded Reload
	add	r12, rcx
	add	r12, r8
	mov	r8, r10
	mov	r10, r14
	mov	rcx, r14
	imul	rcx, rdx
	add	rcx, r12
	imul	r9, qword ptr [rbp - 112]       ## 8-byte Folded Reload
	imul	rdi, rax
	add	rdi, r9
	imul	rbx, rsi
	imul	r11, qword ptr [rbp - 104]      ## 8-byte Folded Reload
	mov	rsi, qword ptr [rbp - 152]      ## 8-byte Reload
	add	r11, rbx
	add	r11, rdi
	imul	rdx, qword ptr [rbp - 72]       ## 8-byte Folded Reload
	add	rdx, r11
	mov	rax, qword ptr [rbp - 144]      ## 8-byte Reload
	mov	rdi, rax
	shr	rdi, 26
	and	eax, 67108863
	mov	rbx, rax
	mov	r9d, edi
	add	r9, r15
	mov	rdi, r9
	shr	rdi, 26
	and	r9d, 67108863
	mov	r12d, edi
	add	r12, r13
	mov	rdi, r12
	shr	rdi, 26
	mov	r13d, edi
	mov	rdi, qword ptr [rbp - 96]       ## 8-byte Reload
	add	r13, rcx
	mov	rcx, r13
	shr	rcx, 26
	mov	r15d, ecx
	add	r15, rdx
	mov	rdx, qword ptr [rbp - 88]       ## 8-byte Reload
	mov	rax, r15
	shr	rax, 26
	lea	r11d, [rax + 4*rax]
	add	r11d, ebx
	mov	r14, qword ptr [rbp - 72]       ## 8-byte Reload
	mov	eax, r11d
	shr	eax, 26
	add	eax, r9d
	and	r12d, 67108863
	and	r13d, 67108863
	and	r15d, 67108863
	and	r11d, 67108863
	add	rsi, 16
	add	rdi, -16
	cmp	rdi, 15
	ja	LBB1_1
LBB1_2:
	test	rdi, rdi
	je	LBB1_9
## %bb.3:
	vxorps	xmm0, xmm0, xmm0
	vmovaps	xmmword ptr [rbp - 64], xmm0
	mov	r9d, edi
	and	r9d, 7
	cmp	rdi, 8
	jae	LBB1_11
## %bb.4:
	xor	ecx, ecx
	mov	rdx, r8
	jmp	LBB1_5
LBB1_11:
	mov	ebx, edi
	and	ebx, 8
	xor	ecx, ecx
	mov	rdx, r8
	.p2align	4, 0x90
LBB1_12:                                ## =>This Inner Loop Header: Depth=1
	movzx	r8d, byte ptr [rsi + rcx]
	mov	byte ptr [rbp + rcx - 64], r8b
	movzx	r8d, byte ptr [rsi + rcx + 1]
	mov	byte ptr [rbp + rcx - 63], r8b
	movzx	r8d, byte ptr [rsi + rcx + 2]
	mov	byte ptr [rbp + rcx - 62], r8b
	movzx	r8d, byte ptr [rsi + rcx + 3]
	mov	byte ptr [rbp + rcx - 61], r8b
	movzx	r8d, byte ptr [rsi + rcx + 4]
	mov	byte ptr [rbp + rcx - 60], r8b
	movzx	r8d, byte ptr [rsi + rcx + 5]
	mov	byte ptr [rbp + rcx - 59], r8b
	movzx	r8d, byte ptr [rsi + rcx + 6]
	mov	byte ptr [rbp + rcx - 58], r8b
	movzx	r8d, byte ptr [rsi + rcx + 7]
	mov	byte ptr [rbp + rcx - 57], r8b
	add	rcx, 8
	cmp	rbx, rcx
	jne	LBB1_12
LBB1_5:
	mov	qword ptr [rbp - 96], rdx       ## 8-byte Spill
	test	r9, r9
	je	LBB1_8
## %bb.6:
	lea	rbx, [rcx + rbp]
	add	rbx, -64
	add	rsi, rcx
	xor	ecx, ecx
	.p2align	4, 0x90
LBB1_7:                                 ## =>This Inner Loop Header: Depth=1
	movzx	r8d, byte ptr [rsi + rcx]
	mov	byte ptr [rbx + rcx], r8b
	inc	rcx
	cmp	r9, rcx
	jne	LBB1_7
LBB1_8:
	mov	byte ptr [rbp + rdi - 64], 1
	mov	esi, 67108863
	and	esi, dword ptr [rbp - 64]
	add	esi, r11d
	mov	edx, dword ptr [rbp - 61]
	shr	edx, 2
	and	edx, 67108863
	add	edx, eax
	mov	r14d, dword ptr [rbp - 58]
	shr	r14d, 4
	and	r14d, 67108863
	add	r14d, r12d
	mov	edi, dword ptr [rbp - 55]
	shr	edi, 6
	add	edi, r13d
	mov	r11d, dword ptr [rbp - 52]
	shr	r11d, 8
	add	r11d, r15d
	mov	rcx, rsi
	mov	rbx, qword ptr [rbp - 72]       ## 8-byte Reload
	imul	rcx, rbx
	mov	r8, r10
	imul	r8, rdx
	add	r8, rcx
	mov	r15, qword ptr [rbp - 96]       ## 8-byte Reload
	mov	rcx, r15
	imul	rcx, r14
	mov	rax, qword ptr [rbp - 80]       ## 8-byte Reload
	mov	r9, rax
	imul	r9, rdi
	add	r9, rcx
	add	r9, r8
	mov	rcx, qword ptr [rbp - 88]       ## 8-byte Reload
	imul	rcx, r11
	add	rcx, r9
	mov	qword ptr [rbp - 88], rcx       ## 8-byte Spill
	mov	rcx, rsi
	mov	r13, qword ptr [rbp - 104]      ## 8-byte Reload
	imul	rcx, r13
	mov	r8, rdx
	imul	r8, rbx
	add	r8, rcx
	mov	rcx, r10
	imul	rcx, r14
	mov	r9, r15
	imul	r9, rdi
	add	r9, rcx
	add	r9, r8
	imul	rax, r11
	add	rax, r9
	mov	qword ptr [rbp - 80], rax       ## 8-byte Spill
	mov	rcx, rsi
	mov	r12, qword ptr [rbp - 128]      ## 8-byte Reload
	imul	rcx, r12
	mov	r8, rdx
	imul	r8, r13
	add	r8, rcx
	mov	rcx, r14
	imul	rcx, rbx
	mov	r9, r10
	imul	r9, rdi
	add	r9, rcx
	add	r9, r8
	imul	r15, r11
	add	r15, r9
	mov	rax, qword ptr [rbp - 112]      ## 8-byte Reload
	imul	rax, rsi
	mov	r9, qword ptr [rbp - 120]       ## 8-byte Reload
	imul	rsi, r9
	mov	rcx, rdx
	imul	rcx, r12
	add	rcx, rsi
	mov	rsi, r14
	imul	rsi, r13
	mov	r8, rdi
	imul	r8, rbx
	add	r8, rsi
	add	r8, rcx
	imul	r10, r11
	add	r10, r8
	imul	rdx, r9
	add	rdx, rax
	imul	r14, r12
	imul	rdi, r13
	add	rdi, r14
	add	rdi, rdx
	imul	r11, rbx
	add	r11, rdi
	mov	rax, qword ptr [rbp - 88]       ## 8-byte Reload
	mov	rcx, rax
	shr	rcx, 26
	and	eax, 67108863
	mov	rsi, rax
	mov	edx, ecx
	add	rdx, qword ptr [rbp - 80]       ## 8-byte Folded Reload
	mov	rcx, rdx
	shr	rcx, 26
	mov	r12d, ecx
	add	r12, r15
	mov	rcx, r12
	shr	rcx, 26
	mov	r13d, ecx
	add	r13, r10
	mov	rcx, r13
	shr	rcx, 26
	mov	r15d, ecx
	add	r15, r11
	mov	rax, r15
	shr	rax, 26
	lea	eax, [rax + 4*rax]
	add	eax, esi
	mov	r11d, eax
	and	edx, 67108863
	shr	eax, 26
	add	eax, edx
	and	r11d, 67108863
	mov	byte ptr [rbp - 64], 0
	mov	byte ptr [rbp - 63], 0
	mov	byte ptr [rbp - 62], 0
	mov	byte ptr [rbp - 61], 0
	mov	byte ptr [rbp - 60], 0
	mov	byte ptr [rbp - 59], 0
	mov	byte ptr [rbp - 58], 0
	mov	byte ptr [rbp - 57], 0
	mov	byte ptr [rbp - 56], 0
	mov	byte ptr [rbp - 55], 0
	mov	byte ptr [rbp - 54], 0
	mov	byte ptr [rbp - 53], 0
	mov	byte ptr [rbp - 52], 0
	mov	byte ptr [rbp - 51], 0
	mov	byte ptr [rbp - 50], 0
	mov	byte ptr [rbp - 49], 0
	and	r12d, 67108863
	and	r13d, 67108863
	and	r15d, 67108863
LBB1_9:
	mov	rcx, qword ptr [rbp - 136]      ## 8-byte Reload
	mov	dword ptr [rcx + 72], r11d
	mov	dword ptr [rcx + 76], eax
	mov	dword ptr [rcx + 80], r12d
	mov	dword ptr [rcx + 84], r13d
	mov	dword ptr [rcx + 88], r15d
	mov	rax, qword ptr [rip + ___stack_chk_guard@GOTPCREL]
	mov	rax, qword ptr [rax]
	cmp	rax, qword ptr [rbp - 48]
	jne	LBB1_13
## %bb.10:
	add	rsp, 120
	pop	rbx
	pop	r12
	pop	r13
	pop	r14
	pop	r15
	pop	rbp
	ret
LBB1_13:
	call	___stack_chk_fail
	.cfi_endproc
                                        ## -- End function
	.globl	_cfx_poly1305_mac               ## -- Begin function cfx_poly1305_mac
	.p2align	4, 0x90
_cfx_poly1305_mac:                      ## @cfx_poly1305_mac
	.cfi_startproc
## %bb.0:
	push	rbp
	.cfi_def_cfa_offset 16
	.cfi_offset rbp, -16
	mov	rbp, rsp
	.cfi_def_cfa_register rbp
	push	r15
	push	r14
	push	r13
	push	r12
	push	rbx
	sub	rsp, 440
	.cfi_offset rbx, -56
	.cfi_offset r12, -48
	.cfi_offset r13, -40
	.cfi_offset r14, -32
	.cfi_offset r15, -24
	mov	qword ptr [rbp - 480], rcx      ## 8-byte Spill
	mov	rax, qword ptr [rip + ___stack_chk_guard@GOTPCREL]
	mov	rax, qword ptr [rax]
	mov	qword ptr [rbp - 48], rax
	mov	eax, dword ptr [rdi]
	mov	ecx, eax
	and	ecx, 67108863
	mov	byte ptr [rbp - 204], al
	mov	byte ptr [rbp - 296], ah
	shr	eax, 16
	mov	byte ptr [rbp - 250], al
	mov	qword ptr [rbp - 400], rcx      ## 8-byte Spill
	mov	eax, ecx
	shr	eax, 24
	mov	byte ptr [rbp - 295], al
	mov	eax, dword ptr [rdi + 3]
	mov	r9d, eax
	shr	r9d, 2
	and	r9d, 67108611
	mov	byte ptr [rbp - 200], r9b
	mov	ecx, eax
	shr	ecx, 10
	mov	byte ptr [rbp - 294], cl
	shr	eax, 18
	mov	byte ptr [rbp - 248], al
	mov	eax, r9d
	shr	eax, 24
	mov	byte ptr [rbp - 293], al
	mov	eax, dword ptr [rdi + 6]
	mov	ecx, eax
	shr	ecx, 4
	mov	byte ptr [rbp - 196], cl
                                        ## kill: def $ecx killed $ecx def $rcx
	and	ecx, 67092735
	mov	byte ptr [rbp - 292], ch
	shr	eax, 20
	mov	byte ptr [rbp - 246], al
	mov	eax, ecx
	mov	r10, rcx
	shr	eax, 24
	mov	byte ptr [rbp - 291], al
	mov	eax, dword ptr [rdi + 9]
	mov	ecx, eax
	shr	ecx, 6
	mov	byte ptr [rbp - 192], cl
	mov	ebx, ecx
	and	ebx, 66076671
	mov	byte ptr [rbp - 290], bh
	mov	ecx, ebx
	shr	ecx, 16
	mov	byte ptr [rbp - 244], cl
	shr	eax, 30
	mov	byte ptr [rbp - 289], al
	mov	eax, dword ptr [rdi + 12]
	mov	ecx, eax
	shr	ecx, 8
	mov	byte ptr [rbp - 188], cl
	mov	r8d, ecx
	and	r8d, 1048575
	shr	eax, 16
	mov	byte ptr [rbp - 288], al
	mov	eax, r8d
	shr	eax, 16
	mov	byte ptr [rbp - 242], al
	mov	qword ptr [rbp - 304], r9       ## 8-byte Spill
	lea	r9, [r9 + 4*r9]
	mov	eax, r9d
	shr	eax, 16
	mov	byte ptr [rbp - 240], al
	mov	eax, r9d
	shr	eax, 24
	mov	byte ptr [rbp - 285], al
	mov	qword ptr [rbp - 360], r10      ## 8-byte Spill
	lea	r10, [r10 + 4*r10]
	mov	eax, r10d
	shr	eax, 16
	mov	byte ptr [rbp - 236], al
	mov	eax, r10d
	shr	eax, 24
	mov	byte ptr [rbp - 281], al
	mov	qword ptr [rbp - 352], rbx      ## 8-byte Spill
	lea	rcx, [rbx + 4*rbx]
	mov	eax, ecx
	shr	eax, 16
	mov	byte ptr [rbp - 232], al
	mov	eax, ecx
	shr	eax, 24
	mov	byte ptr [rbp - 277], al
	mov	qword ptr [rbp - 376], r8       ## 8-byte Spill
	lea	rbx, [r8 + 4*r8]
	mov	eax, ebx
	shr	eax, 16
	mov	byte ptr [rbp - 228], al
	mov	byte ptr [rbp - 287], 0
	mov	rax, r9
	mov	byte ptr [rbp - 136], al
	mov	qword ptr [rbp - 408], r9       ## 8-byte Spill
	mov	byte ptr [rbp - 286], ah
	mov	byte ptr [rbp - 184], 0
	mov	byte ptr [rbp - 284], 0
	mov	byte ptr [rbp - 238], 0
	mov	byte ptr [rbp - 283], 0
	mov	rax, r10
	mov	byte ptr [rbp - 128], al
	mov	qword ptr [rbp - 392], r10      ## 8-byte Spill
	mov	byte ptr [rbp - 282], ah
	mov	byte ptr [rbp - 180], 0
	mov	byte ptr [rbp - 280], 0
	mov	byte ptr [rbp - 234], 0
	mov	byte ptr [rbp - 279], 0
	mov	byte ptr [rbp - 120], cl
	mov	qword ptr [rbp - 384], rcx      ## 8-byte Spill
	mov	byte ptr [rbp - 278], ch
	mov	byte ptr [rbp - 176], 0
	mov	byte ptr [rbp - 276], 0
	mov	byte ptr [rbp - 230], 0
	mov	byte ptr [rbp - 275], 0
	mov	byte ptr [rbp - 112], bl
	mov	qword ptr [rbp - 328], rbx      ## 8-byte Spill
	mov	byte ptr [rbp - 274], bh
	mov	byte ptr [rbp - 273], 0
	mov	byte ptr [rbp - 172], 0
	mov	byte ptr [rbp - 272], 0
	mov	byte ptr [rbp - 226], 0
	mov	byte ptr [rbp - 271], 0
	movzx	eax, byte ptr [rdi + 16]
	movzx	ecx, byte ptr [rdi + 17]
	movzx	r8d, byte ptr [rdi + 18]
	movzx	r9d, byte ptr [rdi + 19]
	mov	r10d, dword ptr [rdi + 16]
	mov	qword ptr [rbp - 448], r10      ## 8-byte Spill
	mov	byte ptr [rbp - 168], al
	mov	byte ptr [rbp - 270], cl
	mov	byte ptr [rbp - 224], r8b
	mov	byte ptr [rbp - 269], r9b
	movzx	eax, byte ptr [rdi + 20]
	movzx	ecx, byte ptr [rdi + 21]
	movzx	r8d, byte ptr [rdi + 22]
	movzx	r9d, byte ptr [rdi + 23]
	mov	r10d, dword ptr [rdi + 20]
	mov	qword ptr [rbp - 456], r10      ## 8-byte Spill
	mov	byte ptr [rbp - 164], al
	mov	byte ptr [rbp - 268], cl
	mov	byte ptr [rbp - 222], r8b
	mov	byte ptr [rbp - 267], r9b
	movzx	eax, byte ptr [rdi + 24]
	movzx	ecx, byte ptr [rdi + 25]
	movzx	r8d, byte ptr [rdi + 26]
	movzx	r9d, byte ptr [rdi + 27]
	mov	r10d, dword ptr [rdi + 24]
	mov	qword ptr [rbp - 464], r10      ## 8-byte Spill
	mov	byte ptr [rbp - 160], al
	mov	byte ptr [rbp - 266], cl
	mov	byte ptr [rbp - 220], r8b
	mov	byte ptr [rbp - 265], r9b
	movzx	eax, byte ptr [rdi + 28]
	movzx	ecx, byte ptr [rdi + 29]
	movzx	r8d, byte ptr [rdi + 30]
	movzx	r9d, byte ptr [rdi + 31]
	mov	edi, dword ptr [rdi + 28]
	mov	qword ptr [rbp - 472], rdi      ## 8-byte Spill
	mov	byte ptr [rbp - 156], al
	mov	byte ptr [rbp - 264], cl
	mov	byte ptr [rbp - 218], r8b
	mov	byte ptr [rbp - 263], r9b
	mov	byte ptr [rbp - 72], 0
	mov	byte ptr [rbp - 77], 0
	mov	byte ptr [rbp - 74], 0
	mov	byte ptr [rbp - 76], 0
	mov	byte ptr [rbp - 84], 0
	mov	byte ptr [rbp - 95], 0
	mov	byte ptr [rbp - 90], 0
	mov	byte ptr [rbp - 94], 0
	mov	byte ptr [rbp - 68], 0
	mov	byte ptr [rbp - 75], 0
	mov	byte ptr [rbp - 88], 0
	mov	byte ptr [rbp - 93], 0
	xor	eax, eax
	cmp	rdx, 16
	jb	LBB2_1
## %bb.2:
	xor	ebx, ebx
	xor	r13d, r13d
	xor	ecx, ecx
	xor	r14d, r14d
	xor	eax, eax
	.p2align	4, 0x90
LBB2_3:                                 ## =>This Inner Loop Header: Depth=1
	mov	qword ptr [rbp - 312], rsi      ## 8-byte Spill
	mov	qword ptr [rbp - 320], rdx      ## 8-byte Spill
	mov	r15d, dword ptr [rsi]
	mov	edx, 67108863
	and	r15d, edx
	add	r15d, ebx
	mov	ebx, dword ptr [rsi + 3]
	mov	r12d, dword ptr [rsi + 6]
	shr	ebx, 2
	and	ebx, 67108863
	add	ebx, r13d
	shr	r12d, 4
	and	r12d, 67108863
	add	r12d, ecx
	mov	r11d, dword ptr [rsi + 9]
	shr	r11d, 6
	add	r11d, r14d
	mov	edx, dword ptr [rsi + 12]
	shr	edx, 8
	lea	r9d, [rax + rdx]
	add	r9d, 16777216
	mov	rax, r15
	mov	rcx, qword ptr [rbp - 400]      ## 8-byte Reload
	imul	rax, rcx
	mov	rsi, rcx
	mov	rdi, qword ptr [rbp - 328]      ## 8-byte Reload
	mov	rdx, rdi
	imul	rdx, rbx
	add	rdx, rax
	mov	r14, qword ptr [rbp - 384]      ## 8-byte Reload
	mov	rax, r14
	imul	rax, r12
	mov	r8, qword ptr [rbp - 392]       ## 8-byte Reload
	mov	r13, r8
	imul	r13, r11
	add	r13, rax
	add	r13, rdx
	mov	r10, qword ptr [rbp - 408]      ## 8-byte Reload
	imul	r10, r9
	add	r10, r13
	mov	rax, r15
	imul	rax, qword ptr [rbp - 304]      ## 8-byte Folded Reload
	mov	rdx, rbx
	imul	rdx, rcx
	add	rdx, rax
	mov	rax, rdi
	imul	rax, r12
	mov	r13, r14
	imul	r13, r11
	add	r13, rax
	add	r13, rdx
	imul	r8, r9
	add	r8, r13
	mov	rdx, r15
	imul	rdx, qword ptr [rbp - 360]      ## 8-byte Folded Reload
	mov	r13, rbx
	mov	rcx, qword ptr [rbp - 304]      ## 8-byte Reload
	imul	r13, rcx
	add	r13, rdx
	mov	rdx, r12
	mov	rax, rsi
	imul	rdx, rsi
	imul	rdi, r11
	add	rdi, rdx
	add	rdi, r13
	imul	r14, r9
	add	r14, rdi
	mov	rdi, r15
	mov	rdx, qword ptr [rbp - 352]      ## 8-byte Reload
	imul	rdi, rdx
	mov	r13, rbx
	mov	rsi, qword ptr [rbp - 360]      ## 8-byte Reload
	imul	r13, rsi
	add	r13, rdi
	mov	rdi, r12
	imul	rdi, rcx
	mov	rcx, r11
	imul	rcx, rax
	add	rcx, rdi
	add	rcx, r13
	mov	r13, qword ptr [rbp - 328]      ## 8-byte Reload
	imul	r13, r9
	add	r13, rcx
	imul	r15, qword ptr [rbp - 376]      ## 8-byte Folded Reload
	imul	rbx, rdx
	add	rbx, r15
	imul	r12, rsi
	imul	r11, qword ptr [rbp - 304]      ## 8-byte Folded Reload
	add	r11, r12
	add	r11, rbx
	imul	r9, rax
	mov	rsi, qword ptr [rbp - 312]      ## 8-byte Reload
	add	r9, r11
	mov	rdx, qword ptr [rbp - 320]      ## 8-byte Reload
	mov	rcx, r10
	shr	rcx, 26
	and	r10d, 67108863
	mov	r15d, ecx
	add	r15, r8
	mov	rax, r15
	shr	rax, 26
	and	r15d, 67108863
	mov	ecx, eax
	add	rcx, r14
	mov	rax, rcx
	shr	rax, 26
	mov	r14d, eax
	add	r14, r13
	mov	rax, r14
	shr	rax, 26
	mov	eax, eax
	add	rax, r9
	mov	rdi, rax
	shr	rdi, 26
	lea	ebx, [rdi + 4*rdi]
	add	ebx, r10d
	mov	r13d, ebx
	shr	r13d, 26
	add	r13d, r15d
	and	ecx, 67108863
	and	r14d, 67108863
	and	eax, 67108863
	and	ebx, 67108863
	add	rsi, 16
	add	rdx, -16
	cmp	rdx, 15
	ja	LBB2_3
	jmp	LBB2_4
LBB2_1:
	xor	r14d, r14d
	xor	ecx, ecx
	xor	r13d, r13d
	xor	ebx, ebx
LBB2_4:
	mov	byte ptr [rbp - 72], bl
	mov	byte ptr [rbp - 77], bh
	mov	edi, ebx
	shr	edi, 16
	mov	byte ptr [rbp - 74], dil
	mov	edi, ebx
	shr	edi, 24
	mov	byte ptr [rbp - 76], dil
	mov	byte ptr [rbp - 84], r13b
	mov	r9d, r13d
	shr	r9d, 8
	mov	byte ptr [rbp - 95], r9b
	mov	r12d, r13d
	shr	r12d, 16
	mov	byte ptr [rbp - 90], r12b
	mov	edi, r13d
	shr	edi, 24
	mov	byte ptr [rbp - 94], dil
	mov	byte ptr [rbp - 68], cl
	mov	r8d, ecx
	shr	r8d, 8
	mov	dword ptr [rbp - 336], r8d      ## 4-byte Spill
	mov	byte ptr [rbp - 75], r8b
	mov	r8d, ecx
	shr	r8d, 16
	mov	dword ptr [rbp - 320], r8d      ## 4-byte Spill
	mov	byte ptr [rbp - 88], r8b
	mov	r8d, ecx
	shr	r8d, 24
	mov	byte ptr [rbp - 93], r8b
	mov	byte ptr [rbp - 96], r14b
	mov	r11d, r14d
	shr	r11d, 8
	mov	byte ptr [rbp - 100], r11b
	mov	r10d, r14d
	shr	r10d, 16
	mov	dword ptr [rbp - 312], r10d     ## 4-byte Spill
	mov	byte ptr [rbp - 98], r10b
	mov	qword ptr [rbp - 368], r14      ## 8-byte Spill
	mov	r10d, r14d
	shr	r10d, 24
	mov	byte ptr [rbp - 99], r10b
	mov	byte ptr [rbp - 80], al
	mov	r14d, eax
	shr	r14d, 8
	mov	byte ptr [rbp - 92], r14b
	mov	r15d, eax
	shr	r15d, 16
	mov	dword ptr [rbp - 344], r15d     ## 4-byte Spill
	mov	byte ptr [rbp - 86], r15b
	mov	r15d, eax
	shr	r15d, 24
	mov	byte ptr [rbp - 91], r15b
	test	rdx, rdx
	je	LBB2_5
## %bb.6:
	vxorps	xmm0, xmm0, xmm0
	vmovaps	xmmword ptr [rbp - 64], xmm0
	mov	eax, edx
	and	eax, 7
	cmp	rdx, 8
	jae	LBB2_14
## %bb.7:
	mov	r8, rdx
	xor	edi, edi
	mov	r14, qword ptr [rbp - 304]      ## 8-byte Reload
	jmp	LBB2_8
LBB2_5:
	mov	rdx, qword ptr [rbp - 368]      ## 8-byte Reload
	mov	rbx, r8
	mov	dword ptr [rbp - 304], r14d     ## 4-byte Spill
	mov	r14d, r11d
	jmp	LBB2_12
LBB2_14:
	mov	r8, rdx
	mov	r9d, edx
	and	r9d, 8
	xor	edi, edi
	mov	r14, qword ptr [rbp - 304]      ## 8-byte Reload
	.p2align	4, 0x90
LBB2_15:                                ## =>This Inner Loop Header: Depth=1
	movzx	r10d, byte ptr [rsi + rdi]
	mov	byte ptr [rbp + rdi - 64], r10b
	movzx	r10d, byte ptr [rsi + rdi + 1]
	mov	byte ptr [rbp + rdi - 63], r10b
	movzx	r10d, byte ptr [rsi + rdi + 2]
	mov	byte ptr [rbp + rdi - 62], r10b
	movzx	r10d, byte ptr [rsi + rdi + 3]
	mov	byte ptr [rbp + rdi - 61], r10b
	movzx	r10d, byte ptr [rsi + rdi + 4]
	mov	byte ptr [rbp + rdi - 60], r10b
	movzx	r10d, byte ptr [rsi + rdi + 5]
	mov	byte ptr [rbp + rdi - 59], r10b
	movzx	r10d, byte ptr [rsi + rdi + 6]
	mov	byte ptr [rbp + rdi - 58], r10b
	movzx	r10d, byte ptr [rsi + rdi + 7]
	mov	byte ptr [rbp + rdi - 57], r10b
	add	rdi, 8
	cmp	r9, rdi
	jne	LBB2_15
LBB2_8:
	test	rax, rax
	mov	rdx, qword ptr [rbp - 368]      ## 8-byte Reload
	je	LBB2_11
## %bb.9:
	lea	r9, [rdi + rbp]
	add	r9, -64
	add	rsi, rdi
	xor	edi, edi
	.p2align	4, 0x90
LBB2_10:                                ## =>This Inner Loop Header: Depth=1
	movzx	r10d, byte ptr [rsi + rdi]
	mov	byte ptr [r9 + rdi], r10b
	inc	rdi
	cmp	rax, rdi
	jne	LBB2_10
LBB2_11:
	mov	byte ptr [rbp + r8 - 64], 1
	mov	eax, 67108863
	and	eax, dword ptr [rbp - 64]
	add	ebx, eax
	mov	eax, dword ptr [rbp - 61]
	shr	eax, 2
	and	eax, 67108863
	add	r13d, eax
	mov	rax, qword ptr [rbp - 328]      ## 8-byte Reload
	mov	rsi, rax
	mov	rdi, rax
	imul	rsi, r13
	mov	qword ptr [rbp - 320], rsi      ## 8-byte Spill
	mov	rsi, rbx
	mov	rax, qword ptr [rbp - 352]      ## 8-byte Reload
	imul	rsi, rax
	mov	qword ptr [rbp - 312], rsi      ## 8-byte Spill
	mov	r10, r13
	mov	qword ptr [rbp - 344], r13      ## 8-byte Spill
	mov	rsi, r13
	imul	r13, rax
	mov	eax, dword ptr [rbp - 58]
	shr	eax, 4
	and	eax, 67108863
	add	ecx, eax
	mov	r15, qword ptr [rbp - 384]      ## 8-byte Reload
	mov	rax, r15
	imul	rax, rcx
	mov	qword ptr [rbp - 416], rax      ## 8-byte Spill
	mov	r8, rdi
	mov	rax, rdi
	imul	rax, rcx
	mov	qword ptr [rbp - 424], rax      ## 8-byte Spill
	mov	rdi, rbx
	mov	rax, qword ptr [rbp - 360]      ## 8-byte Reload
	imul	rdi, rax
	mov	qword ptr [rbp - 440], rdi      ## 8-byte Spill
	imul	rsi, rax
	mov	qword ptr [rbp - 352], rsi      ## 8-byte Spill
	mov	qword ptr [rbp - 432], rcx      ## 8-byte Spill
	mov	qword ptr [rbp - 336], rcx      ## 8-byte Spill
	imul	rcx, rax
	mov	eax, dword ptr [rbp - 55]
	shr	eax, 6
	add	edx, eax
	mov	r9, qword ptr [rbp - 392]       ## 8-byte Reload
	mov	rax, r9
	imul	rax, rdx
	mov	rsi, r15
	imul	rsi, rdx
	mov	rdi, r8
	imul	rdi, rdx
	mov	r12, rdx
	imul	rdx, r14
	add	rdx, rcx
	movzx	ecx, byte ptr [rbp - 91]
	shl	ecx, 24
	movzx	r11d, byte ptr [rbp - 86]
	shl	r11d, 16
	or	r11d, ecx
	movzx	r8d, byte ptr [rbp - 92]
	shl	r8d, 8
	or	r8d, r11d
	movzx	ecx, byte ptr [rbp - 80]
	or	r8d, ecx
	mov	ecx, dword ptr [rbp - 52]
	shr	ecx, 8
	add	r8d, ecx
	mov	rcx, qword ptr [rbp - 376]      ## 8-byte Reload
	imul	rcx, rbx
	add	r13, rcx
	add	rdx, r13
	mov	r14d, dword ptr [rbp - 400]     ## 4-byte Reload
	mov	r11, qword ptr [rbp - 408]      ## 8-byte Reload
	imul	r11, r8
	imul	r9, r8
	imul	r15, r8
	mov	r13, qword ptr [rbp - 328]      ## 8-byte Reload
	imul	r13, r8
	imul	r8, r14
	add	r8, rdx
	mov	rcx, rbx
	imul	rcx, r14
	mov	rdx, qword ptr [rbp - 320]      ## 8-byte Reload
	add	rdx, rcx
	add	rax, qword ptr [rbp - 416]      ## 8-byte Folded Reload
	add	rax, rdx
	add	r11, rax
	mov	rdx, qword ptr [rbp - 304]      ## 8-byte Reload
	imul	rbx, rdx
	imul	r10, r14
	add	r10, rbx
	add	rsi, qword ptr [rbp - 424]      ## 8-byte Folded Reload
	add	rsi, r10
	add	r9, rsi
	mov	rax, r11
	shr	rax, 26
	mov	ecx, eax
	add	rcx, r9
	mov	rax, qword ptr [rbp - 344]      ## 8-byte Reload
	imul	rax, rdx
	mov	rsi, rdx
	add	rax, qword ptr [rbp - 440]      ## 8-byte Folded Reload
	mov	rdx, qword ptr [rbp - 432]      ## 8-byte Reload
	imul	rdx, r14
	add	rdi, rdx
	add	rdi, rax
	add	r15, rdi
	mov	rax, rcx
	shr	rax, 26
	mov	edi, eax
	add	rdi, r15
	mov	rax, qword ptr [rbp - 352]      ## 8-byte Reload
	add	rax, qword ptr [rbp - 312]      ## 8-byte Folded Reload
	mov	rdx, qword ptr [rbp - 336]      ## 8-byte Reload
	imul	rdx, rsi
	imul	r12, r14
	add	r12, rdx
	add	r12, rax
	mov	rsi, r13
	add	rsi, r12
	mov	rax, rdi
	shr	rax, 26
	mov	edx, eax
	add	rdx, rsi
	and	r11d, 67108863
	mov	rax, rdx
	shr	rax, 26
	mov	eax, eax
	add	rax, r8
	mov	rsi, rax
	shr	rsi, 26
	lea	ebx, [rsi + 4*rsi]
	add	ebx, r11d
	mov	r13d, ebx
	mov	byte ptr [rbp - 72], bl
	mov	byte ptr [rbp - 77], bh
	mov	esi, ebx
	shr	ebx, 24
	and	bl, 3
	mov	byte ptr [rbp - 76], bl
	mov	rbx, rdi
	and	ecx, 67108863
	mov	byte ptr [rbp - 68], bl
	mov	byte ptr [rbp - 75], bh
	shr	edi, 16
	mov	dword ptr [rbp - 320], edi      ## 4-byte Spill
	mov	byte ptr [rbp - 88], dil
	shr	ebx, 24
	and	bl, 3
	mov	byte ptr [rbp - 93], bl
	mov	byte ptr [rbp - 96], dl
	mov	r14d, edx
	shr	r14d, 8
	mov	byte ptr [rbp - 100], r14b
	mov	edi, edx
	shr	edi, 16
	mov	dword ptr [rbp - 312], edi      ## 4-byte Spill
	mov	byte ptr [rbp - 98], dil
	mov	r10d, edx
	shr	r10d, 24
	and	r10b, 3
	mov	byte ptr [rbp - 99], r10b
	mov	byte ptr [rbp - 80], al
	mov	edi, eax
	shr	edi, 8
	mov	dword ptr [rbp - 304], edi      ## 4-byte Spill
	mov	byte ptr [rbp - 92], dil
	mov	edi, eax
	shr	edi, 16
	mov	dword ptr [rbp - 344], edi      ## 4-byte Spill
	mov	byte ptr [rbp - 86], dil
	mov	r15d, eax
	shr	r15d, 24
	and	r15b, 3
	mov	byte ptr [rbp - 91], r15b
	shr	r13d, 26
	add	r13d, ecx
	shr	esi, 16
	mov	byte ptr [rbp - 74], sil
	mov	byte ptr [rbp - 84], r13b
	mov	r9d, r13d
	shr	r9d, 8
	mov	byte ptr [rbp - 95], r9b
	mov	r12d, r13d
	shr	r12d, 16
	mov	byte ptr [rbp - 90], r12b
	mov	edi, r13d
	shr	edi, 24
	mov	byte ptr [rbp - 94], dil
	mov	byte ptr [rbp - 64], 0
	mov	byte ptr [rbp - 63], 0
	mov	byte ptr [rbp - 62], 0
	mov	byte ptr [rbp - 61], 0
	mov	byte ptr [rbp - 60], 0
	mov	byte ptr [rbp - 59], 0
	mov	byte ptr [rbp - 58], 0
	mov	byte ptr [rbp - 57], 0
	mov	byte ptr [rbp - 56], 0
	mov	byte ptr [rbp - 55], 0
	mov	byte ptr [rbp - 54], 0
	mov	byte ptr [rbp - 53], 0
	mov	byte ptr [rbp - 52], 0
	mov	byte ptr [rbp - 51], 0
	mov	byte ptr [rbp - 50], 0
	mov	byte ptr [rbp - 49], 0
	movzx	ecx, byte ptr [rbp - 75]
	mov	dword ptr [rbp - 336], ecx      ## 4-byte Spill
	movzx	ecx, byte ptr [rbp - 68]
LBB2_12:
	movzx	esi, bl
	shl	r9d, 8
	movzx	r8d, r9w
	movzx	r9d, r13b
	or	r9d, r8d
	movzx	r8d, r12b
	shl	r8d, 16
	mov	r11d, edi
	and	edi, 3
	shl	edi, 24
	or	edi, r8d
	or	edi, r9d
	shl	esi, 24
	movzx	r8d, byte ptr [rbp - 320]       ## 1-byte Folded Reload
	shl	r8d, 16
	or	r8d, esi
	movzx	esi, byte ptr [rbp - 336]       ## 1-byte Folded Reload
	shl	esi, 8
	or	esi, r8d
	movzx	r8d, cl
	or	r8d, esi
	shr	r11d, 2
	add	r8d, r11d
	movzx	ecx, r10b
	shl	ecx, 24
	movzx	esi, byte ptr [rbp - 312]       ## 1-byte Folded Reload
	shl	esi, 16
	or	esi, ecx
	movzx	ecx, r14b
	shl	ecx, 8
	or	ecx, esi
	movzx	r9d, dl
	or	r9d, ecx
	mov	ecx, r8d
	shr	ecx, 26
	and	r8d, 67108863
	add	r9d, ecx
	mov	ecx, r9d
	shr	ecx, 26
	and	r9d, 67108863
	movzx	edx, r15b
	shl	edx, 24
	movzx	esi, byte ptr [rbp - 344]       ## 1-byte Folded Reload
	shl	esi, 16
	or	esi, edx
	movzx	edx, byte ptr [rbp - 304]       ## 1-byte Folded Reload
	shl	edx, 8
	or	edx, esi
	movzx	esi, al
	or	esi, edx
	add	esi, ecx
	mov	eax, esi
	shr	eax, 26
	lea	eax, [rax + 4*rax]
	movzx	ecx, byte ptr [rbp - 72]
	movzx	edx, byte ptr [rbp - 74]
	movzx	r10d, byte ptr [rbp - 76]
	shl	r10d, 24
	shl	edx, 16
	or	edx, r10d
	movzx	r15d, byte ptr [rbp - 77]
	shl	r15d, 8
	or	r15d, edx
	or	r15d, ecx
	add	r15d, eax
	mov	r14d, r15d
	shr	r14d, 26
	add	r14d, edi
	and	r15d, 67108863
	lea	eax, [r15 + 5]
	mov	ebx, eax
	mov	r12d, eax
	shr	ebx, 26
	add	ebx, r14d
	mov	edx, ebx
	shr	edx, 26
	add	edx, r8d
	mov	r10d, edx
	shr	r10d, 26
	add	r10d, r9d
	mov	eax, r10d
	shr	eax, 26
	mov	edi, esi
	or	edi, -67108864
	add	edi, eax
	mov	r11d, edi
	sar	r11d, 31
	mov	eax, edi
	shr	eax, 31
	dec	eax
	mov	qword ptr [rbp - 304], rax      ## 8-byte Spill
	and	r15d, r11d
	mov	ecx, eax
	and	ecx, 67108863
	mov	eax, r12d
	and	eax, ecx
	or	eax, r15d
	and	r14d, r11d
	and	ebx, ecx
	or	ebx, r14d
	mov	byte ptr [rbp - 72], al
	mov	byte ptr [rbp - 77], ah
	mov	r13d, ebx
	shl	r13d, 26
	or	r13d, eax
	mov	dword ptr [rbp - 328], eax      ## 4-byte Spill
	mov	r15d, eax
	and	r8d, r11d
	and	edx, ecx
	or	edx, r8d
	mov	byte ptr [rbp - 84], bl
	mov	byte ptr [rbp - 95], bh
	mov	r12d, ebx
	mov	r14d, ebx
	shr	ebx, 6
	mov	r8d, edx
	shl	r8d, 20
	or	r8d, ebx
	and	esi, r11d
	and	r11d, r9d
	and	ecx, r10d
	or	ecx, r11d
	mov	byte ptr [rbp - 68], dl
	mov	byte ptr [rbp - 75], dh
	mov	r10d, edx
	mov	r9d, edx
	shr	edx, 12
	mov	ebx, ecx
	shl	ebx, 14
	or	ebx, edx
	mov	rdx, qword ptr [rbp - 304]      ## 8-byte Reload
	and	edx, edi
	and	esi, 67108863
	or	edx, esi
	mov	byte ptr [rbp - 96], cl
	mov	byte ptr [rbp - 100], ch
	mov	byte ptr [rbp - 80], dl
	mov	byte ptr [rbp - 92], dh
	mov	eax, ecx
	mov	esi, ecx
	shr	ecx, 18
	mov	edi, edx
	mov	r11d, edx
	shl	edx, 8
	or	edx, ecx
	shr	r15d, 16
	mov	byte ptr [rbp - 74], r15b
	mov	ecx, dword ptr [rbp - 328]      ## 4-byte Reload
	shr	ecx, 24
	mov	byte ptr [rbp - 76], cl
	shr	r12d, 16
	mov	byte ptr [rbp - 90], r12b
	shr	r14d, 24
	mov	byte ptr [rbp - 94], r14b
	shr	r10d, 16
	mov	byte ptr [rbp - 88], r10b
	shr	r9d, 24
	mov	byte ptr [rbp - 93], r9b
	shr	eax, 16
	mov	byte ptr [rbp - 98], al
	shr	esi, 24
	mov	byte ptr [rbp - 99], sil
	shr	edi, 16
	mov	byte ptr [rbp - 86], dil
	shr	r11d, 24
	mov	byte ptr [rbp - 91], r11b
	add	r13, qword ptr [rbp - 448]      ## 8-byte Folded Reload
	mov	eax, r13d
	shr	eax, 16
	mov	byte ptr [rbp - 216], al
	mov	eax, r13d
	shr	eax, 24
	mov	byte ptr [rbp - 261], al
	add	r8, qword ptr [rbp - 456]       ## 8-byte Folded Reload
	mov	rax, r13
	shr	rax, 32
	add	r8, rax
	mov	eax, r8d
	shr	eax, 16
	mov	byte ptr [rbp - 214], al
	mov	eax, r8d
	shr	eax, 24
	mov	byte ptr [rbp - 259], al
	add	rbx, qword ptr [rbp - 464]      ## 8-byte Folded Reload
	mov	rax, r8
	shr	rax, 32
	add	rbx, rax
	mov	eax, ebx
	shr	eax, 16
	mov	byte ptr [rbp - 212], al
	mov	eax, ebx
	shr	eax, 24
	mov	byte ptr [rbp - 257], al
	add	rdx, qword ptr [rbp - 472]      ## 8-byte Folded Reload
	mov	rax, rbx
	shr	rax, 32
	add	rdx, rax
	mov	rax, rdx
	shr	rax, 32
	mov	byte ptr [rbp - 140], al
	mov	eax, edx
	shr	eax, 8
	mov	byte ptr [rbp - 254], al
	mov	byte ptr [rbp - 256], al
	mov	eax, edx
	shr	eax, 16
	mov	byte ptr [rbp - 208], al
	mov	byte ptr [rbp - 210], al
	mov	eax, edx
	shr	eax, 24
	mov	byte ptr [rbp - 253], al
	mov	byte ptr [rbp - 255], al
	mov	rax, r13
	mov	byte ptr [rbp - 64], al
	mov	byte ptr [rbp - 262], ah
	mov	rcx, r13
	mov	rax, r8
	mov	byte ptr [rbp - 152], al
	mov	byte ptr [rbp - 260], ah
	mov	rsi, r8
	mov	byte ptr [rbp - 148], bl
	mov	byte ptr [rbp - 258], bh
	mov	byte ptr [rbp - 104], dl
	mov	byte ptr [rbp - 252], 0
	mov	byte ptr [rbp - 206], 0
	mov	byte ptr [rbp - 251], 0
	mov	byte ptr [rbp - 144], dl
	mov	rax, qword ptr [rbp - 480]      ## 8-byte Reload
	mov	dword ptr [rax], ecx
	mov	dword ptr [rax + 4], esi
	mov	dword ptr [rax + 8], ebx
	mov	dword ptr [rax + 12], edx
	mov	byte ptr [rbp - 204], 0
	mov	byte ptr [rbp - 296], 0
	mov	byte ptr [rbp - 250], 0
	mov	byte ptr [rbp - 295], 0
	mov	byte ptr [rbp - 200], 0
	mov	byte ptr [rbp - 294], 0
	mov	byte ptr [rbp - 248], 0
	mov	byte ptr [rbp - 293], 0
	mov	byte ptr [rbp - 196], 0
	mov	byte ptr [rbp - 292], 0
	mov	byte ptr [rbp - 246], 0
	mov	byte ptr [rbp - 291], 0
	mov	byte ptr [rbp - 192], 0
	mov	byte ptr [rbp - 290], 0
	mov	byte ptr [rbp - 244], 0
	mov	byte ptr [rbp - 289], 0
	mov	byte ptr [rbp - 188], 0
	mov	byte ptr [rbp - 288], 0
	mov	byte ptr [rbp - 242], 0
	mov	byte ptr [rbp - 287], 0
	mov	byte ptr [rbp - 136], 0
	mov	byte ptr [rbp - 286], 0
	mov	byte ptr [rbp - 240], 0
	mov	byte ptr [rbp - 285], 0
	mov	byte ptr [rbp - 184], 0
	mov	byte ptr [rbp - 284], 0
	mov	byte ptr [rbp - 238], 0
	mov	byte ptr [rbp - 283], 0
	mov	byte ptr [rbp - 128], 0
	mov	byte ptr [rbp - 282], 0
	mov	byte ptr [rbp - 236], 0
	mov	byte ptr [rbp - 281], 0
	mov	byte ptr [rbp - 180], 0
	mov	byte ptr [rbp - 280], 0
	mov	byte ptr [rbp - 234], 0
	mov	byte ptr [rbp - 279], 0
	mov	byte ptr [rbp - 120], 0
	mov	byte ptr [rbp - 278], 0
	mov	byte ptr [rbp - 232], 0
	mov	byte ptr [rbp - 277], 0
	mov	byte ptr [rbp - 176], 0
	mov	byte ptr [rbp - 276], 0
	mov	byte ptr [rbp - 230], 0
	mov	byte ptr [rbp - 275], 0
	mov	byte ptr [rbp - 112], 0
	mov	byte ptr [rbp - 274], 0
	mov	byte ptr [rbp - 228], 0
	mov	byte ptr [rbp - 273], 0
	mov	byte ptr [rbp - 172], 0
	mov	byte ptr [rbp - 272], 0
	mov	byte ptr [rbp - 226], 0
	mov	byte ptr [rbp - 271], 0
	mov	byte ptr [rbp - 72], 0
	mov	byte ptr [rbp - 77], 0
	mov	byte ptr [rbp - 74], 0
	mov	byte ptr [rbp - 76], 0
	mov	byte ptr [rbp - 84], 0
	mov	byte ptr [rbp - 95], 0
	mov	byte ptr [rbp - 90], 0
	mov	byte ptr [rbp - 94], 0
	mov	byte ptr [rbp - 68], 0
	mov	byte ptr [rbp - 75], 0
	mov	byte ptr [rbp - 88], 0
	mov	byte ptr [rbp - 93], 0
	mov	byte ptr [rbp - 96], 0
	mov	byte ptr [rbp - 100], 0
	mov	byte ptr [rbp - 98], 0
	mov	byte ptr [rbp - 99], 0
	mov	byte ptr [rbp - 80], 0
	mov	byte ptr [rbp - 92], 0
	mov	byte ptr [rbp - 86], 0
	mov	byte ptr [rbp - 91], 0
	mov	byte ptr [rbp - 168], 0
	mov	byte ptr [rbp - 270], 0
	mov	byte ptr [rbp - 224], 0
	mov	byte ptr [rbp - 269], 0
	mov	byte ptr [rbp - 164], 0
	mov	byte ptr [rbp - 268], 0
	mov	byte ptr [rbp - 222], 0
	mov	byte ptr [rbp - 267], 0
	mov	byte ptr [rbp - 160], 0
	mov	byte ptr [rbp - 266], 0
	mov	byte ptr [rbp - 220], 0
	mov	byte ptr [rbp - 265], 0
	mov	byte ptr [rbp - 156], 0
	mov	byte ptr [rbp - 264], 0
	mov	byte ptr [rbp - 218], 0
	mov	byte ptr [rbp - 263], 0
	mov	byte ptr [rbp - 64], 0
	mov	byte ptr [rbp - 262], 0
	mov	byte ptr [rbp - 216], 0
	mov	byte ptr [rbp - 261], 0
	mov	byte ptr [rbp - 152], 0
	mov	byte ptr [rbp - 260], 0
	mov	byte ptr [rbp - 214], 0
	mov	byte ptr [rbp - 259], 0
	mov	byte ptr [rbp - 148], 0
	mov	byte ptr [rbp - 258], 0
	mov	byte ptr [rbp - 212], 0
	mov	byte ptr [rbp - 257], 0
	mov	byte ptr [rbp - 144], 0
	mov	byte ptr [rbp - 256], 0
	mov	byte ptr [rbp - 210], 0
	mov	byte ptr [rbp - 255], 0
	mov	byte ptr [rbp - 104], 0
	mov	byte ptr [rbp - 254], 0
	mov	byte ptr [rbp - 208], 0
	mov	byte ptr [rbp - 253], 0
	mov	byte ptr [rbp - 140], 0
	mov	byte ptr [rbp - 252], 0
	mov	byte ptr [rbp - 206], 0
	mov	byte ptr [rbp - 251], 0
	mov	rax, qword ptr [rip + ___stack_chk_guard@GOTPCREL]
	mov	rax, qword ptr [rax]
	cmp	rax, qword ptr [rbp - 48]
	jne	LBB2_16
## %bb.13:
	add	rsp, 440
	pop	rbx
	pop	r12
	pop	r13
	pop	r14
	pop	r15
	pop	rbp
	ret
LBB2_16:
	call	___stack_chk_fail
	.cfi_endproc
                                        ## -- End function
.subsections_via_symbols
