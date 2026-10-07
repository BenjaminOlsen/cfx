BGE provides password-based file encryption through the cfx_bge library
(static or shared), the `cfx bge` command, and the standalone `cfx_bge` tool.
It uses Argon2id and XChaCha20-Poly1305, with a fresh random salt and nonce
for each encryption. The output stores both values for decryption.

Usage (prompts for a password):
  cfx bge -e -i plain.txt -o secret.bge
  cfx bge -d -i secret.bge -o plain.txt

Input/output default to stdin/stdout. Add -a for Base64-armored encryption;
decryption detects armor automatically.

Code layout:
  bge.h               Public library interface.
  cfx_bge.c           Library implementation.
  cfx_bge_internal.h  Private constants and declarations.
  cfx_bge_cli.c       Argument parsing, password prompts, and CLI file I/O.
  cfx_bge_main.c      Standalone executable entry point.

Link consumers against the CMake target cfx::bge and include <bge.h>.

