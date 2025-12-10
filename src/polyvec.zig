//! Vectors of polynomials for Dilithium.
//! PolyVecL has L polynomials, PolyVecK has K polynomials.

const std = @import("std");
const params = @import("params.zig");
const poly_mod = @import("poly.zig");

const K = params.K;
const L = params.L;
const SEEDBYTES = params.SEEDBYTES;
const CRHBYTES = params.CRHBYTES;
const POLYW1_PACKEDBYTES = params.POLYW1_PACKEDBYTES;

const Poly = poly_mod.Poly;

/// Vector of L polynomials.
pub const PolyVecL = struct {
    vec: [L]Poly = [_]Poly{.{}} ** L,

    const Self = @This();

    /// Sample vector with uniformly random coefficients in [-ETA, ETA].
    pub fn uniformEta(self: *Self, seed: *const [CRHBYTES]u8, nonce: u16) void {
        // Validate nonce to prevent overflow
        if (@as(u32, nonce) > @as(u32, std.math.maxInt(u16)) - L) return;
        inline for (&self.vec, 0..) |*p, i| p.uniformEta(seed, nonce + @as(u16, @intCast(i)));
    }

    /// Sample vector with uniformly random coefficients in [-(GAMMA1-1), GAMMA1].
    pub fn uniformGamma1(self: *Self, seed: *const [CRHBYTES]u8, nonce: u16) void {
        // Validate nonce to prevent overflow in multiplication
        if (@as(u32, nonce) > @as(u32, std.math.maxInt(u16)) / L) return;
        inline for (&self.vec, 0..) |*p, i| p.uniformGamma1(seed, @as(u16, L) * nonce + @as(u16, @intCast(i)));
    }

    /// Reduce all coefficients to representatives in [-6283008, 6283008].
    pub fn reduce(self: *Self) void {
        inline for (&self.vec) |*p| p.reduce();
    }

    /// Add vectors: self = u + v. No modular reduction.
    pub fn add(self: *Self, u: *const Self, v: *const Self) void {
        inline for (&self.vec, u.vec, v.vec) |*w, up, vp| w.add(&up, &vp);
    }

    /// Forward NTT of all polynomials in vector.
    /// Output coefficients can be up to 16*Q larger than input.
    pub fn ntt(self: *Self) void {
        inline for (&self.vec) |*p| p.ntt();
    }

    /// Inverse NTT and multiplication by 2^{32} of all polynomials.
    pub fn invnttTomont(self: *Self) void {
        inline for (&self.vec) |*p| p.invnttTomont();
    }

    /// Pointwise multiply all polynomials by a single polynomial with Montgomery reduction.
    /// r[i] = a * v[i] * 2^{-32}
    pub fn pointwisePolyMontgomery(self: *Self, a: *const Poly, v: *const Self) void {
        inline for (&self.vec, v.vec) |*r, vp| r.pointwiseMontgomery(a, &vp);
    }

    /// Pointwise multiply vectors and accumulate with Montgomery reduction.
    /// w = sum(u[i] * v[i]) * 2^{-32}
    pub fn pointwiseAccMontgomery(w: *Poly, u: *const Self, v: *const Self) void {
        w.pointwiseMontgomery(&u.vec[0], &v.vec[0]);

        var t: Poly = .{};
        for (1..L) |i| {
            t.pointwiseMontgomery(&u.vec[i], &v.vec[i]);
            w.add(w, &t);
        }
    }

    /// Check infinity norm of all polynomials in vector.
    /// Returns true if norm of all polynomials is strictly smaller than bound.
    pub fn chknorm(self: *const Self, bound: i32) bool {
        for (self.vec) |p| {
            if (!p.chknorm(bound)) return false;
        }
        return true;
    }
};

/// Vector of K polynomials.
pub const PolyVecK = struct {
    vec: [K]Poly = [_]Poly{.{}} ** K,

    const Self = @This();

    /// Sample vector with uniformly random coefficients in [-ETA, ETA].
    pub fn uniformEta(self: *Self, seed: *const [CRHBYTES]u8, nonce: u16) void {
        // Validate nonce to prevent overflow
        if (@as(u32, nonce) > @as(u32, std.math.maxInt(u16)) - K) return;
        for (&self.vec, 0..) |*p, i| {
            p.uniformEta(seed, nonce + @as(u16, @intCast(i)));
        }
    }

    /// Reduce all coefficients to representatives in [-6283008, 6283008].
    pub fn reduce(self: *Self) void {
        for (&self.vec) |*p| {
            p.reduce();
        }
    }

    /// For all coefficients, add Q if coefficient is negative.
    pub fn caddQ(self: *Self) void {
        for (&self.vec) |*p| {
            p.caddQ();
        }
    }

    /// Add vectors: self = u + v. No modular reduction.
    pub fn add(self: *Self, u: *const Self, v: *const Self) void {
        for (&self.vec, u.vec, v.vec) |*w, up, vp| {
            w.add(&up, &vp);
        }
    }

    /// Subtract vectors: self = u - v. No modular reduction.
    pub fn sub(self: *Self, u: *const Self, v: *const Self) void {
        for (&self.vec, u.vec, v.vec) |*w, up, vp| {
            w.sub(&up, &vp);
        }
    }

    /// Multiply all polynomials by 2^D without modular reduction.
    /// Assumes input coefficients less than 2^{31-D}.
    pub fn shiftl(self: *Self) void {
        for (&self.vec) |*p| {
            p.shiftl();
        }
    }

    /// Forward NTT of all polynomials in vector.
    /// Output coefficients can be up to 16*Q larger than input.
    pub fn ntt(self: *Self) void {
        for (&self.vec) |*p| {
            p.ntt();
        }
    }

    /// Inverse NTT and multiplication by 2^{32} of all polynomials.
    /// Input coefficients need to be less than 2*Q.
    pub fn invnttTomont(self: *Self) void {
        for (&self.vec) |*p| {
            p.invnttTomont();
        }
    }

    /// Pointwise multiply all polynomials by a single polynomial with Montgomery reduction.
    /// r[i] = a * v[i] * 2^{-32}
    pub fn pointwisePolyMontgomery(self: *Self, a: *const Poly, v: *const Self) void {
        for (&self.vec, v.vec) |*r, vp| {
            r.pointwiseMontgomery(a, &vp);
        }
    }

    /// Check infinity norm of all polynomials in vector.
    /// Returns true if norm of all polynomials is strictly smaller than bound.
    pub fn chknorm(self: *const Self, bound: i32) bool {
        for (self.vec) |p| {
            if (!p.chknorm(bound)) return false;
        }
        return true;
    }

    /// For all coefficients a, compute a0, a1 such that a mod^+ Q = a1*2^D + a0
    /// with -2^{D-1} < a0 <= 2^{D-1}.
    pub fn power2Round(v1: *Self, v0: *Self, v: *const Self) void {
        for (&v1.vec, &v0.vec, v.vec) |*p1, *p0, p| {
            Poly.power2Round(p1, p0, &p);
        }
    }

    /// For all coefficients a, compute high and low bits a0, a1 such that
    /// a mod^+ Q = a1*ALPHA + a0.
    pub fn decompose(v1: *Self, v0: *Self, v: *const Self) void {
        for (&v1.vec, &v0.vec, v.vec) |*p1, *p0, p| {
            Poly.polyDecompose(p1, p0, &p);
        }
    }

    /// Compute hint vector. Returns number of 1 bits.
    pub fn makeHint(h: *Self, v0: *const Self, v1: *const Self) u32 {
        var s: u32 = 0;
        for (&h.vec, v0.vec, v1.vec) |*hp, p0, p1| {
            s += Poly.polyMakeHint(hp, &p0, &p1);
        }
        return s;
    }

    /// Use hint vector to correct the high bits of input vector.
    pub fn useHint(w: *Self, u: *const Self, h: *const Self) void {
        for (&w.vec, u.vec, h.vec) |*wp, up, hp| {
            Poly.polyUseHint(wp, &up, &hp);
        }
    }

    /// Pack w1 vector into byte array.
    pub fn packW1(self: *const Self, r: *[@as(usize, K) * POLYW1_PACKEDBYTES]u8) void {
        for (self.vec, 0..) |p, i| {
            p.w1Pack(r[i * POLYW1_PACKEDBYTES ..][0..POLYW1_PACKEDBYTES]);
        }
    }
};

/// Matrix of K x L polynomials (K rows, L columns).
/// Used for the public matrix A in Dilithium.
pub const Mat = struct {
    rows: [K]PolyVecL = [_]PolyVecL{.{}} ** K,

    const Self = @This();

    /// Expand matrix A from seed rho using SHAKE128.
    /// Generates uniformly random coefficients via rejection sampling.
    pub fn expand(self: *Self, rho: *const [SEEDBYTES]u8) void {
        var idx: usize = 0;
        const total = K * L;

        // Process 4 polynomials at a time with parallel SHAKE instances
        while (idx + 4 <= total) : (idx += 4) {
            var nonces: [4]u16 = undefined;
            var polys: [4]*Poly = undefined;

            for (0..4) |k| {
                const flat_idx = idx + k;
                const i = flat_idx / L;
                const j = flat_idx % L;
                nonces[k] = (@as(u16, @intCast(i)) << 8) | @as(u16, @intCast(j));
                polys[k] = &self.rows[i].vec[j];
            }

            uniformX4(polys, rho, nonces);
        }

        // Handle remaining polynomials (0-3)
        while (idx < total) : (idx += 1) {
            const i = idx / L;
            const j = idx % L;
            const nonce: u16 = (@as(u16, @intCast(i)) << 8) | @as(u16, @intCast(j));
            self.rows[i].vec[j].uniform(rho, nonce);
        }
    }

    /// Matrix-vector multiplication with Montgomery reduction.
    /// t = A * v (in NTT domain)
    pub fn pointwiseMontgomery(t: *PolyVecK, self: *const Self, v: *const PolyVecL) void {
        for (&t.vec, self.rows) |*tp, row| {
            PolyVecL.pointwiseAccMontgomery(tp, &row, v);
        }
    }
};

/// Process 4 polynomials in parallel using interleaved SHAKE128 states.
fn uniformX4(
    polys: [4]*Poly,
    rho: *const [SEEDBYTES]u8,
    nonces: [4]u16,
) void {
    // Initialize 4 SHAKE128 states
    var states: [4]std.crypto.hash.sha3.Shake128 = undefined;
    for (0..4) |k| {
        states[k] = std.crypto.hash.sha3.Shake128.init(.{});
        states[k].update(rho);
        states[k].update(std.mem.asBytes(&nonces[k]));
    }

    // Larger buffer reduces squeeze overhead - 840 bytes = 280 coefficients worth
    const BUFLEN = 840;
    var bufs: [4][BUFLEN]u8 = undefined;
    var positions: [4]usize = .{ BUFLEN, BUFLEN, BUFLEN, BUFLEN };
    var coeff_counts: [4]usize = .{ 0, 0, 0, 0 };

    // Track which polynomials are done
    var done: u4 = 0;

    while (done != 0xF) {
        // Squeeze new bytes for states that need them
        for (0..4) |k| {
            if ((done >> @intCast(k)) & 1 == 0 and positions[k] >= BUFLEN - 2) {
                states[k].squeeze(&bufs[k]);
                positions[k] = 0;
            }
        }

        // Process bytes in batches using SIMD extraction
        for (0..4) |k| {
            if ((done >> @intCast(k)) & 1 == 1) continue;

            const remaining_coeffs = params.N - coeff_counts[k];
            const remaining_bytes = BUFLEN - positions[k];
            // Each coefficient needs 3 bytes
            const batch_size = @min(remaining_bytes / 3, remaining_coeffs);

            const accepted = extractCoeffsBatch(
                bufs[k][positions[k]..],
                polys[k].coeffs[coeff_counts[k]..],
                batch_size,
            );

            coeff_counts[k] += accepted;
            positions[k] += batch_size * 3;

            if (coeff_counts[k] >= params.N) {
                done |= @as(u4, 1) << @intCast(k);
            }
        }
    }
}

/// Extract coefficients using optimized SIMD operations.
/// Returns number of accepted coefficients.
fn extractCoeffsBatch(
    src: []const u8,
    dst: []i32,
    max_candidates: usize,
) usize {
    const Q: u32 = 8380417;
    var accepted: usize = 0;
    var i: usize = 0;

    // Process 8 candidates at a time with SIMD
    while (i + 8 <= max_candidates and accepted + 8 <= dst.len) {
        const coeffs = extractEightCoeffs(src[i * 3 ..]);
        const mask = simdReject(coeffs, Q);

        // Use SIMD compress to pack accepted values
        accepted += simdCompress(coeffs, mask, dst[accepted..]);
        i += 8;
    }

    // Scalar fallback for remainder
    while (i < max_candidates and accepted < dst.len) {
        const b0: u32 = src[i * 3 + 0];
        const b1: u32 = src[i * 3 + 1];
        const b2: u32 = src[i * 3 + 2];

        const coeff = b0 | (b1 << 8) | ((b2 & 0x7F) << 16);

        // Constant-time conditional store
        const accept = @intFromBool(coeff < Q);
        dst[accepted] = @intCast(coeff);
        accepted += accept;
        i += 1;
    }

    return accepted;
}

/// Extract 8 23-bit coefficients using SIMD.
inline fn extractEightCoeffs(src: []const u8) @Vector(8, u32) {
    // Load 24 bytes (8 * 3) and extract 23-bit values
    // Use gather pattern optimized for the 3-byte stride

    var result: @Vector(8, u32) = @splat(0);

    // Unrolled extraction - compiler will optimize memory access patterns
    comptime var j: usize = 0;
    inline while (j < 8) : (j += 1) {
        const base = j * 3;
        const b0: u32 = src[base + 0];
        const b1: u32 = src[base + 1];
        const b2: u32 = src[base + 2] & 0x7F;
        result[j] = b0 | (b1 << 8) | (b2 << 16);
    }

    return result;
}

/// SIMD comparison returning mask of values < Q.
inline fn simdReject(coeffs: @Vector(8, u32), q: u32) @Vector(8, bool) {
    const q_vec: @Vector(8, u32) = @splat(q);
    return coeffs < q_vec;
}

/// Compress accepted values and return count.
inline fn simdCompress(
    coeffs: @Vector(8, u32),
    mask: @Vector(8, bool),
    dst: []i32,
) usize {
    var count: usize = 0;

    // Branchless store using mask
    comptime var j: usize = 0;
    inline while (j < 8) : (j += 1) {
        const accept = @intFromBool(mask[j]);
        if (count < dst.len) {
            dst[count] = @intCast(coeffs[j]);
        }
        count += accept;
    }

    return count;
}

// ==================== Tests ====================

test "polyvecl basic" {
    var u = PolyVecL{};
    var v = PolyVecL{};
    var w = PolyVecL{};

    for (&u.vec) |*p| {
        for (&p.coeffs) |*c| {
            c.* = 1;
        }
    }
    for (&v.vec) |*p| {
        for (&p.coeffs) |*c| {
            c.* = 2;
        }
    }

    w.add(&u, &v);
    try std.testing.expectEqual(@as(i32, 3), w.vec[0].coeffs[0]);
}

test "polyveck basic" {
    var u = PolyVecK{};
    var v = PolyVecK{};
    var w = PolyVecK{};

    for (&u.vec) |*p| {
        for (&p.coeffs) |*c| {
            c.* = 5;
        }
    }
    for (&v.vec) |*p| {
        for (&p.coeffs) |*c| {
            c.* = 3;
        }
    }

    w.sub(&u, &v);
    try std.testing.expectEqual(@as(i32, 2), w.vec[0].coeffs[0]);
}

test "polyveck chknorm" {
    var v = PolyVecK{};
    for (&v.vec) |*p| {
        for (&p.coeffs) |*c| {
            c.* = 100;
        }
    }

    try std.testing.expect(v.chknorm(200));
    try std.testing.expect(!v.chknorm(50));
}

test "matrix dimensions" {
    const mat = Mat{};
    try std.testing.expectEqual(@as(usize, K), mat.rows.len);
    try std.testing.expectEqual(@as(usize, L), mat.rows[0].vec.len);
}
