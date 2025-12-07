//! SIMD vector types and utilities for AVX2 acceleration.
//! Provides 256-bit vector operations mapping to AVX2 instructions.

const std = @import("std");
const builtin = @import("builtin");

/// 8-lane i32 vector (256-bit, maps to YMM registers on AVX2).
pub const I32x8 = @Vector(8, i32);

/// 8-lane i64 vector for intermediate products.
pub const I64x8 = @Vector(8, i64);

/// Check if AVX2 is available at compile time for the target.
pub const has_avx2 = blk: {
    if (builtin.cpu.arch != .x86_64 and builtin.cpu.arch != .x86) {
        break :blk false;
    }
    break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.x86.Feature.avx2));
};

/// Check if WASM SIMD is available at compile time.
pub const has_wasm_simd = blk: {
    if (builtin.cpu.arch != .wasm32 and builtin.cpu.arch != .wasm64) {
        break :blk false;
    }
    break :blk builtin.cpu.features.isEnabled(@intFromEnum(std.Target.wasm.Feature.simd128));
};

/// Generic SIMD availability check.
pub const has_simd = has_avx2 or has_wasm_simd;

/// Broadcast a scalar to all lanes.
pub inline fn broadcast(x: i32) I32x8 {
    return @splat(x);
}

/// Load 8 i32s from a slice.
pub inline fn load(ptr: *const [8]i32) I32x8 {
    return ptr.*;
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
