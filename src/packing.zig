//! Packing and unpacking routines for Dilithium keys and signatures.

const std = @import("std");
const params = @import("params.zig");
const poly_mod = @import("poly.zig");
const polyvec_mod = @import("polyvec.zig");

const K = params.K;
const L = params.L;
const N = params.N;
const OMEGA = params.OMEGA;
const SEEDBYTES = params.SEEDBYTES;
const TRBYTES = params.TRBYTES;
const CTILDEBYTES = params.CTILDEBYTES;
const POLYT1_PACKEDBYTES = params.POLYT1_PACKEDBYTES;
const POLYT0_PACKEDBYTES = params.POLYT0_PACKEDBYTES;
const POLYETA_PACKEDBYTES = params.POLYETA_PACKEDBYTES;
const POLYZ_PACKEDBYTES = params.POLYZ_PACKEDBYTES;
const CRYPTO_PUBLICKEYBYTES = params.CRYPTO_PUBLICKEYBYTES;
const CRYPTO_SECRETKEYBYTES = params.CRYPTO_SECRETKEYBYTES;
const CRYPTO_BYTES = params.CRYPTO_BYTES;

const Poly = poly_mod.Poly;
const PolyVecL = polyvec_mod.PolyVecL;
const PolyVecK = polyvec_mod.PolyVecK;

/// Packed public key: pk = (rho, t1).
pub const PublicKey = struct {
    rho: [SEEDBYTES]u8 = [_]u8{0} ** SEEDBYTES,
    t1: PolyVecK = .{},

    const Self = @This();

    /// Pack public key into byte array.
    pub fn pack(self: *const Self, pk: *[CRYPTO_PUBLICKEYBYTES]u8) void {
        @memcpy(pk[0..SEEDBYTES], &self.rho);

        for (self.t1.vec, 0..) |p, i| {
            p.t1Pack(pk[SEEDBYTES + i * POLYT1_PACKEDBYTES ..][0..POLYT1_PACKEDBYTES]);
        }
    }

    /// Unpack public key from byte array.
    pub fn unpack(self: *Self, pk: *const [CRYPTO_PUBLICKEYBYTES]u8) void {
        @memcpy(&self.rho, pk[0..SEEDBYTES]);

        for (&self.t1.vec, 0..) |*p, i| {
            p.t1Unpack(pk[SEEDBYTES + i * POLYT1_PACKEDBYTES ..][0..POLYT1_PACKEDBYTES]);
        }
    }

    /// Create packed public key bytes directly.
    pub fn toBytes(self: *const Self) [CRYPTO_PUBLICKEYBYTES]u8 {
        var pk: [CRYPTO_PUBLICKEYBYTES]u8 = undefined;
        self.pack(&pk);
        return pk;
    }

    /// Create from packed bytes.
    pub fn fromBytes(pk: *const [CRYPTO_PUBLICKEYBYTES]u8) Self {
        var self = Self{};
        self.unpack(pk);
        return self;
    }
};

/// Packed secret key: sk = (rho, key, tr, s1, s2, t0).
pub const SecretKey = struct {
    rho: [SEEDBYTES]u8 = [_]u8{0} ** SEEDBYTES,
    key: [SEEDBYTES]u8 = [_]u8{0} ** SEEDBYTES,
    tr: [TRBYTES]u8 = [_]u8{0} ** TRBYTES,
    s1: PolyVecL = .{},
    s2: PolyVecK = .{},
    t0: PolyVecK = .{},

    const Self = @This();

    // Offsets for packing
    const KEY_OFFSET: usize = SEEDBYTES;
    const TR_OFFSET: usize = KEY_OFFSET + SEEDBYTES;
    const S1_OFFSET: usize = TR_OFFSET + TRBYTES;
    const S2_OFFSET: usize = S1_OFFSET + @as(usize, L) * POLYETA_PACKEDBYTES;
    const T0_OFFSET: usize = S2_OFFSET + @as(usize, K) * POLYETA_PACKEDBYTES;

    /// Pack secret key into byte array.
    pub fn pack(self: *const Self, sk: *[CRYPTO_SECRETKEYBYTES]u8) void {
        @memcpy(sk[0..SEEDBYTES], &self.rho);
        @memcpy(sk[KEY_OFFSET..][0..SEEDBYTES], &self.key);
        @memcpy(sk[TR_OFFSET..][0..TRBYTES], &self.tr);

        for (self.s1.vec, 0..) |p, i| {
            p.etaPack(sk[S1_OFFSET + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
        }

        for (self.s2.vec, 0..) |p, i| {
            p.etaPack(sk[S2_OFFSET + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
        }

        for (self.t0.vec, 0..) |p, i| {
            p.t0Pack(sk[T0_OFFSET + i * POLYT0_PACKEDBYTES ..][0..POLYT0_PACKEDBYTES]);
        }
    }

    /// Unpack secret key from byte array.
    pub fn unpack(self: *Self, sk: *const [CRYPTO_SECRETKEYBYTES]u8) void {
        @memcpy(&self.rho, sk[0..SEEDBYTES]);
        @memcpy(&self.key, sk[KEY_OFFSET..][0..SEEDBYTES]);
        @memcpy(&self.tr, sk[TR_OFFSET..][0..TRBYTES]);

        for (&self.s1.vec, 0..) |*p, i| {
            p.etaUnpack(sk[S1_OFFSET + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
        }

        for (&self.s2.vec, 0..) |*p, i| {
            p.etaUnpack(sk[S2_OFFSET + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
        }

        for (&self.t0.vec, 0..) |*p, i| {
            p.t0Unpack(sk[T0_OFFSET + i * POLYT0_PACKEDBYTES ..][0..POLYT0_PACKEDBYTES]);
        }
    }

    /// Create packed secret key bytes directly.
    pub fn toBytes(self: *const Self) [CRYPTO_SECRETKEYBYTES]u8 {
        var sk: [CRYPTO_SECRETKEYBYTES]u8 = undefined;
        self.pack(&sk);
        return sk;
    }

    /// Create from packed bytes.
    pub fn fromBytes(sk: *const [CRYPTO_SECRETKEYBYTES]u8) Self {
        var self = Self{};
        self.unpack(sk);
        return self;
    }
};

/// Packed signature: sig = (c, z, h).
pub const Signature = struct {
    c: [CTILDEBYTES]u8 = [_]u8{0} ** CTILDEBYTES,
    z: PolyVecL = .{},
    h: PolyVecK = .{},

    const Self = @This();

    // Offsets for packing
    const Z_OFFSET: usize = CTILDEBYTES;
    const H_OFFSET: usize = Z_OFFSET + @as(usize, L) * POLYZ_PACKEDBYTES;

    /// Pack signature into byte array.
    pub fn pack(self: *const Self, sig: *[CRYPTO_BYTES]u8) void {
        @memcpy(sig[0..CTILDEBYTES], &self.c);

        for (self.z.vec, 0..) |p, i| {
            p.zPack(sig[Z_OFFSET + i * POLYZ_PACKEDBYTES ..][0..POLYZ_PACKEDBYTES]);
        }

        // Encode h (sparse representation)
        @memset(sig[H_OFFSET..][0 .. OMEGA + K], 0);

        var k: usize = 0;
        for (self.h.vec, 0..) |p, i| {
            for (p.coeffs, 0..) |coeff, j| {
                if (coeff != 0) {
                    sig[H_OFFSET + k] = @intCast(j);
                    k += 1;
                }
            }
            sig[H_OFFSET + OMEGA + i] = @intCast(k);
        }
    }

    /// Unpack signature from byte array.
    /// Returns error if signature is malformed.
    pub fn unpack(self: *Self, sig: *const [CRYPTO_BYTES]u8) error{MalformedSignature}!void {
        @memcpy(&self.c, sig[0..CTILDEBYTES]);

        for (&self.z.vec, 0..) |*p, i| {
            p.zUnpack(sig[Z_OFFSET + i * POLYZ_PACKEDBYTES ..][0..POLYZ_PACKEDBYTES]);
        }

        // Decode h
        var k: usize = 0;
        for (&self.h.vec, 0..) |*p, i| {
            @memset(&p.coeffs, 0);

            const end = sig[H_OFFSET + OMEGA + i];
            if (end < k or end > OMEGA) {
                return error.MalformedSignature;
            }

            var j: usize = k;
            while (j < end) : (j += 1) {
                // Coefficients must be strictly increasing for strong unforgeability
                if (j > k and sig[H_OFFSET + j] <= sig[H_OFFSET + j - 1]) {
                    return error.MalformedSignature;
                }
                p.coeffs[sig[H_OFFSET + j]] = 1;
            }

            k = end;
        }

        // Extra indices must be zero for strong unforgeability
        for (k..OMEGA) |j| {
            if (sig[H_OFFSET + j] != 0) {
                return error.MalformedSignature;
            }
        }
    }

    /// Create packed signature bytes directly.
    pub fn toBytes(self: *const Self) [CRYPTO_BYTES]u8 {
        var sig: [CRYPTO_BYTES]u8 = undefined;
        self.pack(&sig);
        return sig;
    }

    /// Create from packed bytes.
    /// Returns error if signature is malformed.
    pub fn fromBytes(sig: *const [CRYPTO_BYTES]u8) error{MalformedSignature}!Self {
        var self = Self{};
        try self.unpack(sig);
        return self;
    }
};

// ==================== Standalone Pack Functions (C-compatible API) ====================

/// Pack public key pk = (rho, t1).
pub fn packPk(
    pk: *[CRYPTO_PUBLICKEYBYTES]u8,
    rho: *const [SEEDBYTES]u8,
    t1: *const PolyVecK,
) void {
    @memcpy(pk[0..SEEDBYTES], rho);

    for (t1.vec, 0..) |p, i| {
        p.t1Pack(pk[SEEDBYTES + i * POLYT1_PACKEDBYTES ..][0..POLYT1_PACKEDBYTES]);
    }
}

/// Unpack public key pk = (rho, t1).
pub fn unpackPk(
    rho: *[SEEDBYTES]u8,
    t1: *PolyVecK,
    pk: *const [CRYPTO_PUBLICKEYBYTES]u8,
) void {
    @memcpy(rho, pk[0..SEEDBYTES]);

    for (&t1.vec, 0..) |*p, i| {
        p.t1Unpack(pk[SEEDBYTES + i * POLYT1_PACKEDBYTES ..][0..POLYT1_PACKEDBYTES]);
    }
}

/// Pack secret key sk = (rho, key, tr, s1, s2, t0).
pub fn packSk(
    sk: *[CRYPTO_SECRETKEYBYTES]u8,
    rho: *const [SEEDBYTES]u8,
    tr: *const [TRBYTES]u8,
    key: *const [SEEDBYTES]u8,
    t0: *const PolyVecK,
    s1: *const PolyVecL,
    s2: *const PolyVecK,
) void {
    var offset: usize = 0;

    @memcpy(sk[offset..][0..SEEDBYTES], rho);
    offset += SEEDBYTES;

    @memcpy(sk[offset..][0..SEEDBYTES], key);
    offset += SEEDBYTES;

    @memcpy(sk[offset..][0..TRBYTES], tr);
    offset += TRBYTES;

    for (s1.vec, 0..) |p, i| {
        p.etaPack(sk[offset + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
    }
    offset += @as(usize, L) * POLYETA_PACKEDBYTES;

    for (s2.vec, 0..) |p, i| {
        p.etaPack(sk[offset + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
    }
    offset += @as(usize, K) * POLYETA_PACKEDBYTES;

    for (t0.vec, 0..) |p, i| {
        p.t0Pack(sk[offset + i * POLYT0_PACKEDBYTES ..][0..POLYT0_PACKEDBYTES]);
    }
}

/// Unpack secret key sk = (rho, key, tr, s1, s2, t0).
pub fn unpackSk(
    rho: *[SEEDBYTES]u8,
    tr: *[TRBYTES]u8,
    key: *[SEEDBYTES]u8,
    t0: *PolyVecK,
    s1: *PolyVecL,
    s2: *PolyVecK,
    sk: *const [CRYPTO_SECRETKEYBYTES]u8,
) void {
    var offset: usize = 0;

    @memcpy(rho, sk[offset..][0..SEEDBYTES]);
    offset += SEEDBYTES;

    @memcpy(key, sk[offset..][0..SEEDBYTES]);
    offset += SEEDBYTES;

    @memcpy(tr, sk[offset..][0..TRBYTES]);
    offset += TRBYTES;

    for (&s1.vec, 0..) |*p, i| {
        p.etaUnpack(sk[offset + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
    }
    offset += @as(usize, L) * POLYETA_PACKEDBYTES;

    for (&s2.vec, 0..) |*p, i| {
        p.etaUnpack(sk[offset + i * POLYETA_PACKEDBYTES ..][0..POLYETA_PACKEDBYTES]);
    }
    offset += @as(usize, K) * POLYETA_PACKEDBYTES;

    for (&t0.vec, 0..) |*p, i| {
        p.t0Unpack(sk[offset + i * POLYT0_PACKEDBYTES ..][0..POLYT0_PACKEDBYTES]);
    }
}

/// Pack signature sig = (c, z, h).
pub fn packSig(
    sig: *[CRYPTO_BYTES]u8,
    c: *const [CTILDEBYTES]u8,
    z: *const PolyVecL,
    h: *const PolyVecK,
) void {
    var offset: usize = 0;

    @memcpy(sig[offset..][0..CTILDEBYTES], c);
    offset += CTILDEBYTES;

    for (z.vec, 0..) |p, i| {
        p.zPack(sig[offset + i * POLYZ_PACKEDBYTES ..][0..POLYZ_PACKEDBYTES]);
    }
    offset += L * POLYZ_PACKEDBYTES;

    // Encode h
    @memset(sig[offset..][0 .. OMEGA + K], 0);

    var k: usize = 0;
    for (h.vec, 0..) |p, i| {
        for (p.coeffs, 0..) |coeff, j| {
            if (coeff != 0) {
                sig[offset + k] = @intCast(j);
                k += 1;
            }
        }
        sig[offset + OMEGA + i] = @intCast(k);
    }
}

/// Unpack signature sig = (c, z, h).
/// Returns true on success, false if malformed.
pub fn unpackSig(
    c: *[CTILDEBYTES]u8,
    z: *PolyVecL,
    h: *PolyVecK,
    sig: *const [CRYPTO_BYTES]u8,
) bool {
    var offset: usize = 0;

    @memcpy(c, sig[offset..][0..CTILDEBYTES]);
    offset += CTILDEBYTES;

    for (&z.vec, 0..) |*p, i| {
        p.zUnpack(sig[offset + i * POLYZ_PACKEDBYTES ..][0..POLYZ_PACKEDBYTES]);
    }
    offset += L * POLYZ_PACKEDBYTES;

    // Decode h
    var k: usize = 0;
    for (&h.vec, 0..) |*p, i| {
        @memset(&p.coeffs, 0);

        const end = sig[offset + OMEGA + i];
        if (end < k or end > OMEGA) {
            return false;
        }

        var j: usize = k;
        while (j < end) : (j += 1) {
            // Coefficients must be strictly increasing for strong unforgeability
            if (j > k and sig[offset + j] <= sig[offset + j - 1]) {
                return false;
            }
            p.coeffs[sig[offset + j]] = 1;
        }

        k = end;
    }

    // Extra indices must be zero for strong unforgeability
    for (k..OMEGA) |j| {
        if (sig[offset + j] != 0) {
            return false;
        }
    }

    return true;
}

// ==================== Tests ====================

test "public key pack/unpack roundtrip" {
    var pk1 = PublicKey{};
    for (&pk1.rho, 0..) |*b, i| {
        b.* = @intCast(i);
    }
    for (&pk1.t1.vec[0].coeffs, 0..) |*c, i| {
        c.* = @intCast(i & 0x3FF); // 10-bit values for t1
    }

    const packd = pk1.toBytes();
    const pk2 = PublicKey.fromBytes(&packd);

    try std.testing.expectEqualSlices(u8, &pk1.rho, &pk2.rho);
    try std.testing.expectEqualSlices(i32, &pk1.t1.vec[0].coeffs, &pk2.t1.vec[0].coeffs);
}

test "secret key pack/unpack roundtrip" {
    var sk1 = SecretKey{};
    for (&sk1.rho, 0..) |*b, i| {
        b.* = @intCast(i);
    }
    for (&sk1.key, 0..) |*b, i| {
        b.* = @intCast(i ^ 0xAA);
    }
    for (&sk1.tr, 0..) |*b, i| {
        b.* = @intCast(i ^ 0x55);
    }

    const packd = sk1.toBytes();
    const sk2 = SecretKey.fromBytes(&packd);

    try std.testing.expectEqualSlices(u8, &sk1.rho, &sk2.rho);
    try std.testing.expectEqualSlices(u8, &sk1.key, &sk2.key);
    try std.testing.expectEqualSlices(u8, &sk1.tr, &sk2.tr);
}

test "signature pack/unpack roundtrip" {
    var sig1 = Signature{};
    for (&sig1.c, 0..) |*b, i| {
        b.* = @intCast(i);
    }

    // Set a few hint bits
    sig1.h.vec[0].coeffs[5] = 1;
    sig1.h.vec[0].coeffs[10] = 1;

    const packd = sig1.toBytes();
    const sig2 = try Signature.fromBytes(&packd);

    try std.testing.expectEqualSlices(u8, &sig1.c, &sig2.c);
    try std.testing.expectEqual(@as(i32, 1), sig2.h.vec[0].coeffs[5]);
    try std.testing.expectEqual(@as(i32, 1), sig2.h.vec[0].coeffs[10]);
    try std.testing.expectEqual(@as(i32, 0), sig2.h.vec[0].coeffs[0]);
}

test "signature malformed detection" {
    var bad_sig: [CRYPTO_BYTES]u8 = [_]u8{0} ** CRYPTO_BYTES;

    // Set invalid hint encoding (end marker < previous)
    const h_offset = CTILDEBYTES + L * POLYZ_PACKEDBYTES;
    bad_sig[h_offset + OMEGA] = 5; // First polynomial claims 5 hints
    bad_sig[h_offset + OMEGA + 1] = 3; // Second claims only 3 (less than 5)

    var sig = Signature{};
    const result = sig.unpack(&bad_sig);
    try std.testing.expectError(error.MalformedSignature, result);
}
