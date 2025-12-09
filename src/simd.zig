//! SIMD vector types and utilities for AVX2/AVX-512 acceleration.
//! Provides 256-bit vector operations mapping to AVX2 instructions.
//! Provides 512-bit vector operations mapping to AVX-512 instructions.

const std = @import("std");
const builtin = @import("builtin");

/// 8-lane i32 vector (256-bit, maps to YMM registers on AVX2).
pub const I32x8 = @Vector(8, i32);

/// 8-lane i64 vector for intermediate products.
pub const I64x8 = @Vector(8, i64);

/// 16-lane i32 vector (512-bit, maps to ZMM registers on AVX-512).
pub const I32x16 = @Vector(16, i32);

/// 16-lane i64 vector for intermediate products.
pub const I64x16 = @Vector(16, i64);

/// Check if AVX2 is available at compile time for the target.
pub const has_avx2 = blk: {
    if (builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86) {
        break :blk false;
    }
    break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx2));
};

/// Check if AVX-512F is available at compile time for the target.
pub const has_avx512f = blk: {
    if (builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86) {
        break :blk false;
    }
    break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx512f));
};

/// Check if WASM SIMD is available at compile time.
pub const has_wasm_simd = blk: {
    if (builtin.cpu.arch != .wasm32 and builtin.cpu.arch != .wasm64) {
        break :blk false;
    }
    break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.wasm.Feature.simd128));
};

/// Check if ARM Neon is available at compile time.
pub const has_neon = blk: {
    if (builtin.cpu.arch == .aarch64) {
        break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.aarch64.Feature.neon));
    }
    if (builtin.cpu.arch == .arm) {
        break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.arm.Feature.neon));
    }
    break :blk false;
};

/// Check if RISC-V Vector extension is available at compile time.
pub const has_riscv_v = blk: {
    if (builtin.cpu.arch.isRISCV()) {
        break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.riscv.Feature.v));
    }
    break :blk false;
};

/// Generic SIMD availability check.
pub const has_simd = has_avx2 or has_avx512f or has_wasm_simd or has_neon or has_riscv_v;

/// Broadcast a scalar to all lanes.
pub inline fn broadcast(x: i32) I32x8 {
    return @splat(x);
}

/// Load 8 i32s from a slice.
pub inline fn load(ptr: *const [8]i32) I32x8 {
    return ptr.*;
}

pub inline fn loadU(ptr: *const [8]i32) I32x8 {
    return @as(*align(@alignOf(i32)) const I32x8, @ptrCast(ptr)).*;
}

/// Store 8 i32s to a slice.
pub inline fn store(ptr: *[8]i32, v: I32x8) void {
    ptr.* = v;
}

/// Widen i32x8 to i64x8 for multiplication.
pub inline fn widen(v: I32x8) I64x8 {
    return @intCast(v);
}

/// Truncate i64x8 to i32x8 (low 32 bits).
pub inline fn truncate64(v: I64x8) I32x8 {
    return @truncate(v);
}

/// Arithmetic right shift by 32 for i64x8.
pub inline fn shr32(v: I64x8) I32x8 {
    const shifted = v >> @splat(32);
    return @truncate(shifted);
}

/// Element-wise wrapping multiply i32 * i32 = i32 (low 32 bits).
pub inline fn mulLo(a: I32x8, b: I32x8) I32x8 {
    return a *% b;
}

/// Element-wise full multiply i32 * i32 = i64.
pub inline fn mulWide(a: I32x8, b: I32x8) I64x8 {
    return widen(a) * widen(b);
}

/// Compute absolute value of i32x8 (Bitwise trick: (v ^ mask) - mask).
pub inline fn abs(v: I32x8) I32x8 {
    // Arithmetic shift right to create a mask of all 1s (if neg) or 0s (if pos)
    const mask = v >> @splat(31);
    return (v ^ mask) -% mask;
}

/// Check if any element in vector 'v' is greater than or equal to 'scalar'.
/// Maps to VPCMPGTD + VPMOVMSKB or similar on AVX2.
pub inline fn anyGe(v: I32x8, scalar: i32) bool {
    const threshold: I32x8 = @splat(scalar);
    // resulting bool vector
    const pred = v >= threshold;
    // Horizontal OR reduction
    return @reduce(.Or, pred);
}

test "simd broadcast and load" {
    const v = broadcast(42);
    const arr: [8]i32 = v;
    for (arr) |x| {
        try std.testing.expectEqual(@as(i32, 42), x);
    }
}

test "simd arithmetic" {
    const a = broadcast(10);
    const b = broadcast(3);
    const sum: [8]i32 = a + b;
    const diff: [8]i32 = a - b;
    try std.testing.expectEqual(@as(i32, 13), sum[0]);
    try std.testing.expectEqual(@as(i32, 7), diff[0]);
}

test "avx512 broadcast and load" {
    const v = broadcast16(42);
    const arr: [16]i32 = v;
    for (arr) |x| {
        try std.testing.expectEqual(@as(i32, 42), x);
    }
}

test "avx512 arithmetic" {
    const a = broadcast16(10);
    const b = broadcast16(3);
    const sum: [16]i32 = a + b;
    const diff: [16]i32 = a - b;
    try std.testing.expectEqual(@as(i32, 13), sum[0]);
    try std.testing.expectEqual(@as(i32, 7), diff[0]);
}

test "avx512 absolute value" {
    const v: I32x16 = .{ -5, 3, -10, 7, -1, 1, -20, 15, -3, 8, -12, 6, -9, 4, -7, 11 };
    const abs_v = abs16(v);
    const result: [16]i32 = abs_v;
    const expected = [16]i32{ 5, 3, 10, 7, 1, 1, 20, 15, 3, 8, 12, 6, 9, 4, 7, 11 };
    for (result, expected) |r, e| {
        try std.testing.expectEqual(e, r);
    }
}

test "avx512 any greater than or equal" {
    const v: I32x16 = .{ 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75 };
    try std.testing.expect(anyGe16(v, 10)); // Should be true (10 >= 10)
    try std.testing.expect(!anyGe16(v, 100)); // Should be false (no element >= 100)
}

// ------------------------------------------------------------------------
// AVX-512 16-lane operations (512-bit ZMM registers)
// ------------------------------------------------------------------------

/// Broadcast a scalar to all 16 lanes.
pub inline fn broadcast16(x: i32) I32x16 {
    return @splat(x);
}

/// Load 16 i32s from a slice.
pub inline fn load16(ptr: *const [16]i32) I32x16 {
    return ptr.*;
}

pub inline fn loadU16(ptr: *const [16]i32) I32x16 {
    return @as(*align(@alignOf(i32)) const I32x16, @ptrCast(ptr)).*;
}

/// Store 16 i32s to a slice.
pub inline fn store16(ptr: *[16]i32, v: I32x16) void {
    ptr.* = v;
}

/// Widen i32x16 to i64x16 for multiplication.
pub inline fn widen16(v: I32x16) I64x16 {
    return @intCast(v);
}

/// Truncate i64x16 to i32x16 (low 32 bits).
pub inline fn truncate64_16(v: I64x16) I32x16 {
    return @truncate(v);
}

/// Element-wise wrapping multiply i32 * i32 = i32 (low 32 bits).
pub inline fn mulLo16(a: I32x16, b: I32x16) I32x16 {
    return a *% b;
}

/// Element-wise full multiply i32 * i32 = i64.
pub inline fn mulWide16(a: I32x16, b: I32x16) I64x16 {
    return widen16(a) * widen16(b);
}

/// Compute absolute value of i32x16 (Bitwise trick: (v ^ mask) - mask).
pub inline fn abs16(v: I32x16) I32x16 {
    // Arithmetic shift right to create a mask of all 1s (if neg) or 0s (if pos)
    const mask = v >> @splat(31);
    return (v ^ mask) -% mask;
}

/// Check if any element in 16-lane vector 'v' is greater than or equal to 'scalar'.
pub inline fn anyGe16(v: I32x16, scalar: i32) bool {
    const threshold: I32x16 = @splat(scalar);
    // resulting bool vector
    const pred = v >= threshold;
    // Horizontal OR reduction
    return @reduce(.Or, pred);
}
