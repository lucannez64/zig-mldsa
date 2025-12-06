const params = @import("params.zig");
const sign = @import("sign.zig");
pub const DILITHIUM_MODE: params.DilithiumMode = .mode3;
pub const signing_mode: sign.SigningMode = .deterministic;
