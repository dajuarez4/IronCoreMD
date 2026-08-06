# Production Guide

## Recommended workflow

The easiest workflow is **Canva or PowerPoint for animation**, followed by **CapCut or DaVinci Resolve for voice and final export**.

1. Create a 16:9 project at 3840 × 2160 pixels.
2. Build seven scenes using the timings in `storyboard.md`.
3. Record the narration in a quiet room, standing 15–20 cm from the microphone.
4. Remove long pauses and background noise before adding animation.
5. Match scene changes to narration phrases instead of reading to a finished video.
6. Add subtle sound effects only if they do not compete with the voice.
7. If music is used, choose licensed instrumental music and keep it roughly 18–24 dB below the narration.
8. Add burned-in English captions for accessibility.

## Export settings

- Container: MP4
- Video codec: H.264
- Resolution: 3840 × 2160
- Frame rate: 30 fps
- Color: RGB / Rec.709
- Audio: AAC, 48 kHz, 256–320 kb/s
- Duration: ideally 58–60 seconds
- Maximum file size: 5 GB
- Orientation: landscape 16:9

## Recording advice

- Aim for one confident take rather than sounding overly rehearsed.
- Emphasize “accurate but expensive,” “turn each structure into a graph,” and “challenge it with physics.”
- Pause briefly after the first sentence and before the final sentence.
- Do not read the numerical value too quickly: say “about ten millielectronvolts per atom.”

## Assets already available

- Poster preview: `../poster/poster_preview.png`
- Graph-kernel parity plot: `../poster/figures/graph_kernel_parity.png`
- Magnetic/non-magnetic phonons: `../poster/figures/bcc_magnetic_phonons.png`
- BCC free-energy plot: `../poster/figures/bcc_free_energy.png`

The Earth cutaway, lattice cartoons, supercomputer, and graph-transformation scenes should be created as separate wide images with no embedded text. Add labels in the video editor so they remain sharp at 4K.

## Generated master

`SC26_IronCoreMD_cartoon_4K_silent.mp4` is a 59.93-second silent 4K master timed to `narration_timing.srt`. Import the MP4 into your editor, add the recorded narration, adjust individual scene boundaries only if your delivery differs from the timing guide, and export again using the settings above.
