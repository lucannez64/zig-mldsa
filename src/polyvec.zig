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
        for (&self.vec, 0..) |*p, i| {
            p.uniformEta(seed, nonce + @as(u16, @intCast(i)));
        }
    }

    /// Sample vector with uniformly random coefficients in [-(GAMMA1-1), GAMMA1].
    pub fn uniformGamma1(self: *Self, seed: *const [CRHBYTES]u8, nonce: u16) void {
        for (&self.vec, 0..) |*p, i| {
            p.uniformGamma1(seed, @as(u16, L) * nonce + @as(u16, @intCast(i)));
        }
    }

    /// Reduce all coefficients to representatives in [-6283008, 6283008].
    pub fn reduce(self: *Self) void {
        for (&self.vec) |*p| {
            p.reduce();
        }
    }

    /// Add vectors: self = u + v. No modular reduction.
    pub fn add(self: *Self, u: *const Self, v: *const Self) void {
        for (&self.vec, u.vec, v.vec) |*w, up, vp| {
            w.add(&up, &vp);
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

    /// Pointwise multiply vectors and accumulate with Montgomery reduction.
    /// w = sum(u[i] * v[i]) * 2^{-32}
    pub fn pointwiseAccMontgomery(w: *Poly, u: *const Self, v: *const Self) void {
        const simd = @import("simd.zig");
        if (simd.has_simd) {
            const reduce_avx2 = @import("reduce_avx2.zig");
            var i: usize = 0;
            while (i < params.N) : (i += 8) {
                var sum = simd.broadcast(0);
                inline for (0..L) |j| {
                    const va = simd.load(@ptrCast(&u.vec[j].coeffs[i]));
                    const vb = simd.load(@ptrCast(&v.vec[j].coeffs[i]));
                    const prod = reduce_avx2.pointwiseMont8(va, vb);
                    sum = sum + prod;
                }
                simd.store(@ptrCast(&w.coeffs[i]), sum);
            }
        } else {
            w.pointwiseMontgomery(&u.vec[0], &v.vec[0]);

            var t: Poly = .{};
            for (1..L) |i| {
                t.pointwiseMontgomery(&u.vec[i], &v.vec[i]);
                w.add(w, &t);
            }
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
        const total = K * L;
        var idx: usize = 0;

        // Process 4 polynomials at a time using parallel SHAKE
        while (idx + 4 <= total) : (idx += 4) {
            const idx0 = idx / L;
            const j0 = idx % L;
            const idx1 = (idx + 1) / L;
            const j1 = (idx + 1) % L;
            const idx2 = (idx + 2) / L;
            const j2 = (idx + 2) % L;
            const idx3 = (idx + 3) / L;
            const j3 = (idx + 3) % L;

            const nonce0: u16 = (@as(u16, @intCast(idx0)) << 8) + @as(u16, @intCast(j0));
            const nonce1: u16 = (@as(u16, @intCast(idx1)) << 8) + @as(u16, @intCast(j1));
            const nonce2: u16 = (@as(u16, @intCast(idx2)) << 8) + @as(u16, @intCast(j2));
            const nonce3: u16 = (@as(u16, @intCast(idx3)) << 8) + @as(u16, @intCast(j3));

            Poly.uniform_4x(&self.rows[idx0].vec[j0], &self.rows[idx1].vec[j1], &self.rows[idx2].vec[j2], &self.rows[idx3].vec[j3], rho, .{ nonce0, nonce1, nonce2, nonce3 });
        }

        // Handle remaining polynomials
        while (idx < total) : (idx += 1) {
            const i = idx / L;
            const j = idx % L;
            const nonce: u16 = (@as(u16, @intCast(i)) << 8) + @as(u16, @intCast(j));
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
