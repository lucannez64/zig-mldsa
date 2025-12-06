//! Dilithium parameter definitions.
//! See FIPS 204 / Dilithium specification for details.

const config = @import("config.zig");

pub const DilithiumMode = enum(u8) { mode2 = 2, mode3 = 3, mode5 = 5 };

// Change this or import from config to select the security level.
pub const DILITHIUM_MODE: DilithiumMode = config.DILITHIUM_MODE;

// Common parameters
pub const SEEDBYTES = 32;
pub const CRHBYTES = 64;
pub const TRBYTES = 64;
pub const RNDBYTES = 32;
pub const N = 256;
pub const Q: i32 = 8380417;
pub const D = 13;
pub const ROOT_OF_UNITY = 1753;

// Mode-dependent parameters
pub const K: u8 = switch (DILITHIUM_MODE) {
    .mode2 => 4,
    .mode3 => 6,
    .mode5 => 8,
};

pub const L: u8 = switch (DILITHIUM_MODE) {
    .mode2 => 4,
    .mode3 => 5,
    .mode5 => 7,
};

pub const ETA: u8 = switch (DILITHIUM_MODE) {
    .mode2 => 2,
    .mode3 => 4,
    .mode5 => 2,
};

pub const TAU: u8 = switch (DILITHIUM_MODE) {
    .mode2 => 39,
    .mode3 => 49,
    .mode5 => 60,
};

pub const BETA: u16 = switch (DILITHIUM_MODE) {
    .mode2 => 78,
    .mode3 => 196,
    .mode5 => 120,
};

pub const GAMMA1: i32 = switch (DILITHIUM_MODE) {
    .mode2 => 1 << 17,
    .mode3, .mode5 => 1 << 19,
};

pub const GAMMA2: i32 = switch (DILITHIUM_MODE) {
    .mode2 => (Q - 1) / 88,
    .mode3, .mode5 => (Q - 1) / 32,
};

pub const OMEGA: u8 = switch (DILITHIUM_MODE) {
    .mode2 => 80,
    .mode3 => 55,
    .mode5 => 75,
};

pub const CTILDEBYTES: u8 = switch (DILITHIUM_MODE) {
    .mode2 => 32,
    .mode3 => 48,
    .mode5 => 64,
};

// Packing sizes
pub const POLYT1_PACKEDBYTES = 320;
pub const POLYT0_PACKEDBYTES = 416;
pub const POLYVECH_PACKEDBYTES = OMEGA + K;

pub const POLYZ_PACKEDBYTES: u16 = if (GAMMA1 == (1 << 17)) 576 else 640;
pub const POLYW1_PACKEDBYTES: u8 = if (GAMMA2 == (Q - 1) / 88) 192 else 128;
pub const POLYETA_PACKEDBYTES: u8 = if (ETA == 2) 96 else 128;

// Crypto sizes
pub const CRYPTO_PUBLICKEYBYTES = SEEDBYTES + @as(u32, K) * POLYT1_PACKEDBYTES;

pub const CRYPTO_SECRETKEYBYTES = 2 * SEEDBYTES +
    TRBYTES +
    @as(u32, L) * POLYETA_PACKEDBYTES +
    @as(u32, K) * POLYETA_PACKEDBYTES +
    @as(u32, K) * POLYT0_PACKEDBYTES;

pub const CRYPTO_BYTES = CTILDEBYTES +
    @as(u32, L) * POLYZ_PACKEDBYTES +
    POLYVECH_PACKEDBYTES;
