# Engine selection — 2026-10-28

The selected base is **FRAKTALFAKTUR_RH_LAB_GPU_v1_8_6**, the most developed
recovered implementation in the accepted native C++17 / Direct3D11 engine
lineage. This selection is based on recovered source bytes and fresh CPU
execution, rather than version-number ordering across unrelated branches.

| Recovered package | Scope | Role |
| --- | --- | --- |
| RH_LAB_GPU v1.8.0 | 18 registered entries | Frozen ancestor |
| RH_LAB_GPU v1.8.4 | 26 entries | Eight integrated complex-dynamics extensions |
| RH_LAB_GPU v1.8.5 | 29 entries | Three finite QQRH Hilbert views |
| **RH_LAB_GPU v1.8.6** | **32 entries** | Adds projective metric, kernel Newton, kernel feedback and a point inspector |
| FRONTIER_FASTBOOT v1.8.1 | 18 entries | Startup optimization sibling |
| FRONTIER_TRACE v1.9.0 | 23 entries | Fastboot sibling with five transfer-trace drivers |
| NATIVE_CATALOG v0.24.0 | 73 descriptors, 245 saved views | Broader historical catalogue, older infrastructure |
| ELLIPTIC_HILBERT_ATLAS v0.3.0 RH_MASTER | 15 modes | Separate spectral/resolvent CPU viewer |

## Source continuity

The selected archive preserves all files of the v1.8.0 ancestor; 66 of the
ancestor's 76 files are byte-identical. It preserves the 29 preceding entries
from v1.8.5. Its catalogue combines 18 RH candidates with 14 extra drivers.
All 39 source filenames of Native Catalog v0.24.0 occur in this base: 33 are
identical, six are modified, and the selected base adds 32 source files.

Native Catalog remains broader in exposed families. The current catalogue
intentionally filters descriptors to the registered RH and extra IDs; copying
old descriptors into it does not reactivate all historical families. This
import does not claim to merge every sibling's features.

FRONTIER_TRACE v1.9.0 is a Direct3D11 sibling. It is distinct from the
standalone OpenGL RH_LAB_GPU v1.9.x replacement rejected in the project
history. Its unique trace drivers and deferred UI refresh are candidates for
a later deliberate merge. Its launcher should not replace the current one
wholesale because compiler discovery and error handling differ.

## Provenance and preservation

Imported archive: `FRAKTALFAKTUR_RH_LAB_GPU_v1_8_6.zip`

SHA-256: `98db03286d9b72d118e658464169dccdb44442ce99d5766fc494d5ef7e8b83f0`

The 98 archive files, including source, launchers, documentation, previews,
historical receipts and `MANIFEST.sha256`, are preserved byte-for-byte at the
repository root. The manifest covers 97 other archive files and verified
successfully. New repository documentation, Makefile integration, portable tests
and the rendering adapter are additive and are not covered by that historical
manifest. The archive version remains 1.8.6; the adapter is a repository addition.

## Fresh validation

On 2026-10-28, GCC 13.3 / C++17 on Linux passed all 11 portable CPU suites.
The 14 extension rendering hashes exactly reproduced the archived v1.8.6
receipt, with serial/four-worker agreement. The 18 original candidate tests,
classic integration, lab integration and focus integration also passed.
See `validation/portable_cpu_2026-10-28.log` and `tests/cpu_self_tests.cpp`.

The repository Makefile also built both portable executables successfully.
The new CLI adapter passed 137 invocations: all 32 registered renderers,
serial/parallel equality, saved-view round trips, 16 preset renders and 22
invalid-input cases. A 640 by 480 kernel Newton BMP was rendered and visually
inspected. These validate the additive adapter without changing engine source.

The Windows-dependent engine-controls suite, native Windows compilation,
interactive GUI/inspector, HLSL compilation and hardware D3D11 execution have
not been validated by this import. Historical receipts remain historical.
Finite numerical visualizations do not certify RH or full QQRH positivity.
