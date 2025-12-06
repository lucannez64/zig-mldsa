//! Montgomery reduction and related finite field operations.
//! For the Dilithium/Kyber cryptographic schemes.

const params = @import("params.zig");

const Q = params.Q;
// q^(-1) mod 2^32
const QINV: i32 = 58728449;
// 2^32 % Q
const MONT: i32 = -4186625;

/// Montgomery reduction.
///
/// For finite field element `a` with -2^31 * Q <= a <= Q * 2^31,
/// computes r ≡ a * 2^-32 (mod Q) such that -Q < r < Q.
pub fn montgomeryReduce(a: i64) i32 {
    const t_low: i32 = @truncate(a);
    const t: i64 = @as(i64, t_low) *% QINV;
    const t_cast: i32 = @truncate(t);
    const product: i64 = @as(i64, t_cast) *% Q;
    return @intCast((a - product) >> 32);
}

/// Reduction modulo Q.
///
/// For finite field element `a` with a <= 2^31 - 2^22 - 1,
/// computes r ≡ a (mod Q) such that -6283008 <= r <= 6283008.
pub fn reduce32(a: i32) i32 {
    const t = (a + (1 << 22)) >> 23;
    return a - t *% Q;
}

/// Conditional add Q.
///
/// Adds Q if input coefficient is negative.
pub fn caddq(a: i32) i32 {
    return a + ((a >> 31) & Q);
}

/// Freeze to standard representative.
///
/// For finite field element `a`, computes standard representative r = a mod+ Q.
pub fn freeze(a: i32) i32 {
    const reduced = reduce32(a);
    return caddq(reduced);
}
