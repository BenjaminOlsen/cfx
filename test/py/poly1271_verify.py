#!/usr/bin/env python3
"""
poly1271 reference implementation for verification.
slow but obviously correct.
"""

M127 = (1 << 127) - 1

def clamp_r(r_bytes):
    r = bytearray(r_bytes)
    r[3]  &= 0x0F
    r[7]  &= 0x0F
    r[11] &= 0x0F
    r[15] &= 0x07
    r[4]  &= 0xFC
    r[8]  &= 0xFC
    r[12] &= 0xFC
    return bytes(r)

def le_bytes_to_int(b):
    return int.from_bytes(b, 'little')

def int_to_le_bytes(n, length):
    return n.to_bytes(length, 'little')

def poly1271(message, key):
    """key = r[0:16] || s[16:32], returns 16-byte tag"""
    assert len(key) == 32

    r_bytes = clamp_r(key[0:16])
    s_bytes = key[16:32]

    r = le_bytes_to_int(r_bytes)
    s = le_bytes_to_int(s_bytes)
    assert r < (1 << 127)

    h = 0
    for i in range(0, len(message), 15):
        block = message[i:i+15]
        n = le_bytes_to_int(block + b'\x01')
        h = ((h + n) * r) % M127

    tag_int = (h + s) % (1 << 128)
    return int_to_le_bytes(tag_int, 16)

def test_vector(name, key_hex, msg_hex, expected_tag_hex):
    key = bytes.fromhex(key_hex)
    msg = bytes.fromhex(msg_hex) if msg_hex else b''
    expected = bytes.fromhex(expected_tag_hex)
    result = poly1271(msg, key)

    if result == expected:
        print(f"PASS {name}: {result.hex()}")
        return True
    print(f"FAIL {name}")
    print(f"  Expected: {expected.hex()}")
    print(f"  Got:      {result.hex()}")
    return False

def generate_test_vectors():
    import hashlib

    def make_key(seed):
        return hashlib.sha256(seed.encode()).digest()

    vectors = []
    lengths = [0, 1, 7, 14, 15, 16, 29, 30, 31, 44, 45, 46, 100, 255, 256, 1000]

    for length in lengths:
        key = make_key(f"key_for_len_{length}")
        msg = bytes([(i * 7 + 13) & 0xFF for i in range(length)])
        tag = poly1271(msg, key)
        vectors.append({
            'name': f'len_{length}',
            'key': key.hex(),
            'msg': msg.hex(),
            'tag': tag.hex(),
            'msg_len': length
        })

    return vectors

def main():
    print("=" * 60)
    print("poly1271 reference")
    print("=" * 60)
    print()

    all_passed = True

    all_passed &= test_vector(
        "empty",
        "0102030405060708090a0b0c0d0e0f10" + "1112131415161718191a1b1c1d1e1f20",
        "",
        "1112131415161718191a1b1c1d1e1f20"
    )

    all_passed &= test_vector(
        "single_block",
        "000102030405060708090a0b0c0d0e0f" + "101112131415161718191a1b1c1d1e1f",
        "48656c6c6f2c20506f6c7931323700",
        "44e9b6e7abf4d5a45e07f3d83178ea1f"
    )

    msg_100 = ''.join(f'{i:02x}' for i in range(100))
    key_100 = ''.join(f'{(i*3+7)&0xff:02x}' for i in range(32))
    all_passed &= test_vector("multi_block_100", key_100, msg_100, "c79db66d6682ce7f1ea6662144fa29b2")

    print()
    vectors = generate_test_vectors()
    for v in vectors:
        print(f"len {v['msg_len']:4d}: {v['tag']}")

    print()
    print("PASSED" if all_passed else "FAILED")
    return all_passed

if __name__ == '__main__':
    import sys
    sys.exit(0 if main() else 1)
