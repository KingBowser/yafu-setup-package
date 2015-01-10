# mp_limb_t mulredc1_13(mp_limb_t * z, const mp_limb_t x, const mp_limb_t * y,
#                 const mp_limb_t *m, mp_limb_t inv_m);
#
# Linux:   z: %rdi, x: %rsi, y: %rdx, m: %rcx, inv_m: %r8
#          Needs %rbx, %rsp, %rbp, %r12-%r15 restored
# Windows: z: %rcx, x: %rdx, y: %r8,  m: %r9, inv_m: 28(%rsp)
#          Needs %rbx, %rbp, %rdi, %rsi, %r12...%15 restored



include(`config.m4')

ifdef(`WINDOWS64_ABI',
`define(`Y_PARAM', `%r8')dnl
define(`INVM_PARAM',`72(%rsp)')dnl'
,
`define(`Y_PARAM', `%rdx')dnl
define(`INVM_PARAM',`%r8')dnl'
)dnl
	TEXT
.align 64 # Opteron L1 code cache line is 64 bytes long
	GLOBL GSYM_PREFIX`'mulredc1_13
	TYPE(GSYM_PREFIX`'mulredc1_`'13,`function')

# Implements multiplication and REDC for one input numbers of LENGTH words
# and a multiplier of one word
ifdef(`WINDOWS64_ABI', `# Uses Windows ABI', `# Uses Linux ABI')

# Values that are referenced only once in the loop over j go into r8 .. r14,
# In the inner loop (over j), tmp, x[i], y, m, and u are constant.
# tmp[j], tmp[j+1], tmp[j+2] are updated frequently. These 8 values
# stay in registers and are referenced as
# YP = y, MP = m, 
# X = x, T0 = tmp[j], T1 = tmp[j+1], CY = carry

define(`T0', `%rsi')dnl
define(`T1', `%rbx')dnl
define(`CY', `%rcx')dnl
define(`CYl', `%ecx')dnl
define(`CYb', `%cl')dnl
define(`X', `%r14')dnl		# register that holds x value
define(`U', `%r11')dnl
define(`YP', `%r9')dnl		# register that points to the y array
define(`MP', `%r10')dnl		# register that points to the m array
define(`ZP', `%rdi')dnl		# register that holds z

`#' Register vars: `T0' = T0, `T1' = T1, `CY' = CY, `X' = X, `U' = U
`#'                `YP' = YP, `MP' = MP

GSYM_PREFIX`'mulredc1_13:


#########################################################################
# i = 0 pass
#########################################################################

`#' register values at loop entry: YP = y, MP = m

# We need to compute u

	movq	(Y_PARAM), %rax		# rax = y[0] (time critical, do first)
	pushq	%rbx
	pushq	%r14
ifdef(`WINDOWS64_ABI',
`	pushq	%rsi
	pushq	%rdi
	movq	%r9, MP			# store m in MP
	movq    Y_PARAM, YP
	movq	%rcx, ZP
	movq	%rdx, X'
,
`	movq	Y_PARAM, YP
	movq	%rcx, MP
	movq    %rsi, X		# store x in X
	# ZP is same as passed in'
)

	xorl	CYl, CYl		# set %CY to 0

	mulq	X			# rdx:rax = y[0] * x

	movq 	%rax, T0		# Move low word of product to T0
	movq	%rdx, T1		# Move high word of product to T1

	imulq	INVM_PARAM, %rax	# %rax = ((x[i]*y[0]+tmp[0])*invm)%2^64
	movq	%rax, U			# this is the new u value

	mulq	(MP)			# multipy u*m[0]
	addq	%rax, T0		# Now %T0 = 0, need not be stored
	movq	8(YP), %rax		# Fetch y[1]
	adcq	%rdx, T1		# 
	setc	CYb
	# CY:T1:T0 <= 2*(2^64-1)^2 <= 2^2*128 - 4*2^64 + 2, hence
	# CY:T1 <= 2*2^64 - 4

ifdef(`WANT_ASSERT', `
        pushf
	testq	T0, T0
	jz	assert1
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
assert1:
	popf
')
define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 1
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	8(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 0(ZP)	# Store T0 in z[1-1]
	movq	16(YP), %rax	# Fetch y[j+1] = y[2] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 2
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	16(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 8(ZP)	# Store T0 in z[2-1]
	movq	24(YP), %rax	# Fetch y[j+1] = y[3] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 3
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	24(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 16(ZP)	# Store T0 in z[3-1]
	movq	32(YP), %rax	# Fetch y[j+1] = y[4] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 4
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	32(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 24(ZP)	# Store T0 in z[4-1]
	movq	40(YP), %rax	# Fetch y[j+1] = y[5] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 5
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	40(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 32(ZP)	# Store T0 in z[5-1]
	movq	48(YP), %rax	# Fetch y[j+1] = y[6] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 6
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	48(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 40(ZP)	# Store T0 in z[6-1]
	movq	56(YP), %rax	# Fetch y[j+1] = y[7] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 7
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	56(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 48(ZP)	# Store T0 in z[7-1]
	movq	64(YP), %rax	# Fetch y[j+1] = y[8] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 8
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	64(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 56(ZP)	# Store T0 in z[8-1]
	movq	72(YP), %rax	# Fetch y[j+1] = y[9] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 9
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	72(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 64(ZP)	# Store T0 in z[9-1]
	movq	80(YP), %rax	# Fetch y[j+1] = y[10] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 10
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	80(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 72(ZP)	# Store T0 in z[10-1]
	movq	88(YP), %rax	# Fetch y[j+1] = y[11] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 11
# Register values at entry: 
# %rax = y[j], X = x, U = u
# T0 = value to store in tmp[j], T1 undefined 
# CY = carry into T1 (is <= 2)
# We have CY:T1 <= 2 * 2^64 - 2

	movq	CY, T1		# T1 = CY <= 1

	# Here, T1:T0 <= 2*2^64 - 2
	mulq	X		# y[j] * x
	# rdx:rax <= (2^64-1)^2 <= 2^128 - 2*2^64 + 1
	addq	%rax, T0	# Add low word to T0
	movq	88(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	# T1:T0 <= 2^128 - 2*2^64 + 1 + 2*2^64 - 2 <= 2^128 - 1, no carry!
ifdef(`WANT_ASSERT', `
	jnc	1f
	lea     _GLOBAL_OFFSET_TABLE_(%rip), %rbx # if we do PIC code, we 
		# need to set rbx; if not, it doesnt hurt
	call	GSYM_PREFIX`'abort@plt
1:
',`')
	
	mulq	U		# m[j]*u
	# rdx:rax <= 2^128 - 2*2^64 + 1, T1:T0 <= 2^128 - 1
	addq	T0, %rax	# Add T0 and low word
	movq	%rax, 80(ZP)	# Store T0 in z[11-1]
	movq	96(YP), %rax	# Fetch y[j+1] = y[12] into %rax
	adcq	%rdx, T1	# Add high word with carry to T1
	setc	CYb		# CY <= 1
	# CY:T1:T0 <= 2^128 - 1 + 2^128 - 2*2^64 + 1 <=
	#             2 * 2^128 - 2*2^64 ==> CY:T1 <= 2 * 2^64 - 2

define(`TT', defn(`T0'))dnl
define(`T0', defn(`T1'))dnl
define(`T1', defn(`TT'))dnl
undefine(`TT')dnl
`#' Now `T0' = T0, `T1' = T1


# Pass for j = 12. Don't fetch new data from y[j+1].

	movq	CY, T1		# T1 = CY <= 1
	
	mulq	X		# y[j] * x[i]
	addq	%rax, T0	# Add low word to T0
	movq	96(MP), %rax	# Fetch m[j] into %rax
	adcq	%rdx, T1 	# Add high word with carry to T1
	mulq    U		# m[j]*u
	addq	%rax, T0	# Add low word to T0
	movq	T0, 88(ZP)	# Store T0 in z[j-1]
	adcq	%rdx, T1	# Add high word with carry to T1
	movq	T1, 96(ZP)	# Store T1 in tmp[j]
	setc	CYb		# %CY <= 1

	movq	CY, %rax	# use carry as return value
ifdef(`WINDOWS64_ABI',
`	popq	%rdi
	popq	%rsi
') dnl
	popq	%r14
	popq	%rbx
	ret
