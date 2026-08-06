from __future__ import annotations

import math
import os
from pathlib import Path

import imageio.v2 as imageio
import numpy as np
from PIL import Image, ImageDraw, ImageEnhance, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parent
POSTER = ROOT.parent / "poster"
OUTPUT = ROOT / "SC26_IronCoreMD_cartoon_4K_silent.mp4"
WIDTH, HEIGHT = 1920, 1080
FPS = 15
DURATION = 60

ORANGE = (245, 128, 37)
DARK_BROWN = (86, 43, 18)
NAVY = (8, 33, 54)
BLUE = (0, 91, 145)
LIGHT_BLUE = (100, 194, 229)
CREAM = (255, 247, 239)
WHITE = (255, 255, 255)
GRAY = (70, 76, 82)
RED = (156, 31, 67)

FONT_REGULAR = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_ROUNDED = "/System/Library/Fonts/Supplemental/Arial Rounded Bold.ttf"


def font(size: int, bold: bool = False, rounded: bool = False):
    path = FONT_ROUNDED if rounded else FONT_BOLD if bold else FONT_REGULAR
    return ImageFont.truetype(path, size=size)


def ease(value: float) -> float:
    value = max(0.0, min(1.0, value))
    return value * value * (3.0 - 2.0 * value)


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def background(top=NAVY, bottom=(17, 67, 91)) -> Image.Image:
    y = np.linspace(0.0, 1.0, HEIGHT, dtype=np.float32)[:, None, None]
    top_arr = np.array(top, dtype=np.float32)[None, None, :]
    bottom_arr = np.array(bottom, dtype=np.float32)[None, None, :]
    pixels = top_arr * (1.0 - y) + bottom_arr * y
    pixels = np.repeat(pixels, WIDTH, axis=1).astype(np.uint8)
    return Image.fromarray(pixels, "RGB")


def text_center(draw, xy, text, size, fill=WHITE, bold=True, anchor="mm"):
    draw.text(xy, text, font=font(size, bold=bold), fill=fill, anchor=anchor, align="center",
              stroke_width=1 if bold else 0, stroke_fill=fill)


def wrap(draw, text, max_width, text_font):
    words = text.split()
    lines = []
    current = ""
    for word in words:
        trial = word if not current else current + " " + word
        if draw.textbbox((0, 0), trial, font=text_font)[2] <= max_width:
            current = trial
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)
    return "\n".join(lines)


def header(draw, title, subtitle=None):
    draw.rounded_rectangle((55, 42, 1865, 150), radius=24, fill=(255, 255, 255, 235), outline=ORANGE, width=5)
    draw.text((92, 73), title, font=font(42, bold=True), fill=DARK_BROWN, anchor="lm")
    if subtitle:
        draw.text((1825, 96), subtitle, font=font(24, bold=True), fill=BLUE, anchor="rm")


def footer(draw, label, progress):
    draw.rounded_rectangle((70, 968, 1850, 1040), radius=18, fill=(5, 23, 38, 225), outline=ORANGE, width=3)
    draw.text((105, 1003), label, font=font(28, bold=True), fill=WHITE, anchor="lm")
    bar_left, bar_right = 1390, 1810
    draw.rounded_rectangle((bar_left, 992, bar_right, 1014), radius=10, fill=(100, 110, 120))
    draw.rounded_rectangle((bar_left, 992, int(lerp(bar_left, bar_right, progress)), 1014), radius=10, fill=ORANGE)


def draw_atom(draw, x, y, radius, color=ORANGE, outline=WHITE):
    draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=color, outline=outline, width=max(2, radius // 7))
    highlight = max(2, radius // 4)
    draw.ellipse((x - radius * 0.45, y - radius * 0.48, x - radius * 0.45 + highlight, y - radius * 0.48 + highlight), fill=(255, 222, 178))


def earth_scene(progress: float) -> Image.Image:
    image = background((3, 17, 34), (8, 38, 58))
    draw = ImageDraw.Draw(image, "RGBA")
    for index in range(95):
        x = (index * 197 + 83) % WIDTH
        y = (index * 83 + 41) % 860
        radius = 1 + (index % 3)
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=(255, 255, 255, 140 + (index % 4) * 25))

    zoom = ease(progress)
    cx = int(lerp(960, 720, zoom))
    cy = 555
    radius = int(lerp(330, 445, zoom))
    layers = [
        (1.00, (30, 103, 158), (126, 196, 103)),
        (0.83, (198, 101, 55), None),
        (0.59, (232, 151, 55), None),
        (0.31, (255, 205, 75), None),
    ]
    for scale, color, land in layers:
        r = int(radius * scale)
        draw.ellipse((cx - r, cy - r, cx + r, cy + r), fill=color, outline=(255, 255, 255, 180), width=4)
        if land:
            draw.polygon([(cx-r*.7,cy-r*.2),(cx-r*.25,cy-r*.55),(cx+r*.05,cy-r*.3),(cx-r*.1,cy+r*.02)], fill=land)
    core_r = int(radius * 0.30)
    for glow in range(6, 0, -1):
        gr = core_r + glow * 18
        draw.ellipse((cx-gr, cy-gr, cx+gr, cy+gr), outline=(255, 180, 30, 35), width=18)
    text_center(draw, (cx, cy), "Fe", 72, fill=DARK_BROWN)
    draw.line((cx, cy - core_r, 1270, 340), fill=ORANGE, width=6)
    draw.rounded_rectangle((1250, 255, 1795, 430), radius=26, fill=(255, 247, 239, 245), outline=ORANGE, width=5)
    draw.text((1285, 295), "EARTH'S INNER CORE", font=font(34, bold=True), fill=DARK_BROWN)
    draw.text((1285, 348), "Extreme pressure + temperature", font=font(27, bold=True), fill=BLUE)
    draw.text((1285, 390), "Dominant element: iron", font=font(25), fill=GRAY)

    if progress > 0.42:
        alpha = int(255 * ease((progress - 0.42) / 0.35))
        for idx, (name, ox) in enumerate([("BCC", 1275), ("FCC", 1450), ("HCP", 1625)]):
            draw.rounded_rectangle((ox, 520, ox + 145, 700), radius=20, fill=(8, 33, 54, alpha), outline=(245, 128, 37, alpha), width=4)
            draw.text((ox + 72, 668), name, font=font(24, bold=True), fill=(255, 255, 255, alpha), anchor="mm")
            points = [(ox+35,555),(ox+110,555),(ox+35,625),(ox+110,625),(ox+72,590)]
            for px, py in points:
                draw_atom(draw, px, py, 14, color=(245,128,37,alpha), outline=(255,255,255,alpha))
    header(draw, "WHAT HAPPENS 5,000+ KM BENEATH OUR FEET?", "IRON UNDER EXTREME CONDITIONS")
    footer(draw, "Earth's hidden laboratory", progress)
    return image


def hpc_scene(progress: float) -> Image.Image:
    image = background((6, 28, 48), (10, 73, 101))
    draw = ImageDraw.Draw(image, "RGBA")
    header(draw, "HIGH-PERFORMANCE COMPUTING + FIRST PRINCIPLES", "DENSITY FUNCTIONAL THEORY")
    for rack in range(4):
        x = 95 + rack * 230
        y = 235
        draw.rounded_rectangle((x, y, x + 185, 830), radius=18, fill=(20, 29, 39), outline=LIGHT_BLUE, width=5)
        draw.rectangle((x + 25, y + 25, x + 160, y + 90), fill=(35, 50, 64))
        for row in range(7):
            yy = y + 125 + row * 61
            draw.rounded_rectangle((x + 22, yy, x + 163, yy + 40), radius=7, fill=(43, 58, 70), outline=(95, 115, 130), width=2)
            for led in range(5):
                glow = int(130 + 125 * ((math.sin(progress * 15 + rack + row + led) + 1) / 2))
                draw.ellipse((x + 36 + led*24, yy + 13, x + 46 + led*24, yy + 23), fill=(40, glow, 120))
    draw.text((520, 875), "SUPERCOMPUTER", font=font(29, bold=True), fill=WHITE, anchor="mm")

    cube_left, cube_top, cube_size = 1110, 280, 500
    draw.rounded_rectangle((1005, 220, 1775, 860), radius=34, fill=(255, 247, 239, 242), outline=ORANGE, width=6)
    draw.text((1390, 265), "QUANTUM-MECHANICAL IRON", font=font(33, bold=True), fill=DARK_BROWN, anchor="mm")
    pts = []
    for ix in range(4):
        for iy in range(4):
            x = cube_left + ix * cube_size / 3
            y = cube_top + iy * cube_size / 3
            pts.append((x, y))
    for x, y in pts:
        draw_atom(draw, int(x), int(y), 25, color=ORANGE, outline=DARK_BROWN)
    for ix in range(4):
        for iy in range(4):
            if ix < 3:
                draw.line((cube_left+ix*cube_size/3, cube_top+iy*cube_size/3, cube_left+(ix+1)*cube_size/3, cube_top+iy*cube_size/3), fill=(0,78,122,150), width=4)
            if iy < 3:
                draw.line((cube_left+ix*cube_size/3, cube_top+iy*cube_size/3, cube_left+ix*cube_size/3, cube_top+(iy+1)*cube_size/3), fill=(0,78,122,150), width=4)
    pulse = 0.55 + 0.45 * math.sin(progress * math.pi * 6) ** 2
    draw.rounded_rectangle((1130, 775, 1650, 830), radius=18, fill=(0, 78, 122, 220))
    draw.rectangle((1145, 790, int(1145 + 490 * progress), 815), fill=(245, 128, 37, 255))
    draw.text((1390, 900), "Structure • Energy • Stability", font=font(29, bold=True), fill=WHITE, anchor="mm")
    footer(draw, "DFT: quantum-mechanical accuracy", progress)
    return image


def limits_scene(progress: float) -> Image.Image:
    image = background((49, 22, 15), (110, 48, 18))
    draw = ImageDraw.Draw(image, "RGBA")
    header(draw, "THE COMPUTATIONAL BOTTLENECK", "ACCURATE — BUT EXPENSIVE")
    draw.rounded_rectangle((120, 235, 865, 845), radius=38, fill=(255,247,239,245), outline=ORANGE, width=7)
    draw.text((492, 305), "SMALL SYSTEMS", font=font(38, bold=True), fill=DARK_BROWN, anchor="mm")
    for i in range(5):
        for j in range(5):
            draw_atom(draw, 260 + i*115, 420 + j*75, 22, color=ORANGE, outline=DARK_BROWN)
    draw.rounded_rectangle((1040, 235, 1785, 845), radius=38, fill=(255,247,239,245), outline=ORANGE, width=7)
    draw.text((1412, 305), "SHORT TIME SCALES", font=font(38, bold=True), fill=DARK_BROWN, anchor="mm")
    center=(1412,565); radius=190
    draw.ellipse((center[0]-radius,center[1]-radius,center[0]+radius,center[1]+radius), fill=WHITE, outline=BLUE, width=12)
    for tick in range(12):
        angle=2*math.pi*tick/12
        x1=center[0]+math.sin(angle)*(radius-18); y1=center[1]-math.cos(angle)*(radius-18)
        x2=center[0]+math.sin(angle)*(radius-42); y2=center[1]-math.cos(angle)*(radius-42)
        draw.line((x1,y1,x2,y2),fill=GRAY,width=6)
    angle=progress*math.pi*6
    draw.line((center[0],center[1],center[0]+math.sin(angle)*120,center[1]-math.cos(angle)*120),fill=RED,width=12)
    draw.line((center[0],center[1],center[0]+math.sin(angle/12)*85,center[1]-math.cos(angle/12)*85),fill=BLUE,width=14)
    draw.ellipse((center[0]-14,center[1]-14,center[0]+14,center[1]+14),fill=DARK_BROWN)
    draw.text((960, 900), "We need a faster model without losing the physics.", font=font(34, bold=True), fill=WHITE, anchor="mm")
    footer(draw, "The scale gap", progress)
    return image


def ml_scene(progress: float) -> Image.Image:
    image = background((5, 37, 54), (0, 92, 120))
    draw = ImageDraw.Draw(image, "RGBA")
    header(draw, "LEARN FROM FIRST-PRINCIPLES DATA", "MACHINE-LEARNING INTERATOMIC POTENTIAL")
    for card in range(5):
        x=100+card*170; y=300+(card%2)*220
        draw.rounded_rectangle((x,y,x+135,y+160),radius=18,fill=(255,247,239,235),outline=ORANGE,width=4)
        for k in range(6):
            ax=x+28+(k%3)*40; ay=y+43+(k//3)*48
            draw_atom(draw,ax,ay,12,color=ORANGE,outline=DARK_BROWN)
        draw.text((x+68,y+138),f"DFT {card+1}",font=font(18,bold=True),fill=DARK_BROWN,anchor="mm")
    draw.polygon([(1250,465),(1100,380),(1100,430),(950,430),(950,500),(1100,500),(1100,550)],fill=ORANGE)
    draw.rounded_rectangle((1280,245,1790,780),radius=45,fill=(255,247,239,245),outline=ORANGE,width=7)
    draw.text((1535,315),"LEARNED POTENTIAL",font=font(38,bold=True),fill=DARK_BROWN,anchor="mm")
    nodes=[(1390,450),(1535,405),(1670,455),(1450,585),(1620,600)]
    edges=[(0,1),(1,2),(0,3),(1,3),(1,4),(2,4),(3,4)]
    for a,b in edges: draw.line((*nodes[a],*nodes[b]),fill=BLUE,width=8)
    for index,(x,y) in enumerate(nodes):
        draw.ellipse((x-34,y-34,x+34,y+34),fill=ORANGE if index%2==0 else LIGHT_BLUE,outline=DARK_BROWN,width=5)
    draw.text((1535,710),"Fast energy + force predictions",font=font(27,bold=True),fill=BLUE,anchor="mm")
    footer(draw,"DFT data → scalable interatomic model",progress)
    return image


def graph_scene(progress: float) -> Image.Image:
    image = background((10, 30, 44), (18, 67, 82))
    draw = ImageDraw.Draw(image, "RGBA")
    header(draw, "FROM ATOMIC STRUCTURE TO GRAPH SPACE", "GRAPH-KERNEL SIMILARITY")
    left_nodes=[(210,350),(390,300),(565,360),(245,575),(445,540),(620,610),(360,760),(590,785)]
    right_nodes=[(1280,310),(1480,290),(1680,380),(1270,580),(1510,540),(1690,660),(1390,790),(1600,810)]
    for i,(x,y) in enumerate(left_nodes):
        draw_atom(draw,x,y,31,color=ORANGE if i%3 else LIGHT_BLUE,outline=WHITE)
    for a in range(len(left_nodes)):
        for b in range(a+1,len(left_nodes)):
            if math.dist(left_nodes[a],left_nodes[b])<255:
                draw.line((*left_nodes[a],*left_nodes[b]),fill=(255,255,255,95),width=4)
    draw.text((420,890),"ATOMIC ENVIRONMENT",font=font(31,bold=True),fill=WHITE,anchor="mm")
    arrow_x=int(lerp(810,1040,ease(progress)))
    draw.polygon([(760,500),(980,500),(980,430),(1120,550),(980,670),(980,600),(760,600)],fill=ORANGE)
    for a in range(len(right_nodes)):
        for b in range(a+1,len(right_nodes)):
            if math.dist(right_nodes[a],right_nodes[b])<300:
                draw.line((*right_nodes[a],*right_nodes[b]),fill=LIGHT_BLUE,width=7)
    for i,(x,y) in enumerate(right_nodes):
        draw.ellipse((x-36,y-36,x+36,y+36),fill=ORANGE if i%3 else LIGHT_BLUE,outline=WHITE,width=5)
    draw.text((1490,890),"LABELED GRAPH",font=font(31,bold=True),fill=WHITE,anchor="mm")
    draw.rounded_rectangle((735,760,1135,850),radius=22,fill=(255,247,239,240),outline=ORANGE,width=4)
    draw.text((935,805),"Kᵢⱼ = k(Gᵢ, Gⱼ)",font=font(34,bold=True),fill=DARK_BROWN,anchor="mm")
    footer(draw,"Graph kernels compare local atomic environments",progress)
    return image


def fit_image(source: Image.Image, box, zoom=1.0):
    x0,y0,x1,y1=box; bw=x1-x0; bh=y1-y0
    scale=max(bw/source.width,bh/source.height)*zoom
    resized=source.resize((int(source.width*scale),int(source.height*scale)),Image.Resampling.LANCZOS)
    left=max(0,(resized.width-bw)//2); top=max(0,(resized.height-bh)//2)
    return resized.crop((left,top,left+bw,top+bh))


def result_scene(progress: float) -> Image.Image:
    image=Image.new("RGB",(WIDTH,HEIGHT),CREAM)
    draw=ImageDraw.Draw(image,"RGBA")
    source=Image.open(POSTER/"figures"/"graph_kernel_parity.png").convert("RGB")
    plot=fit_image(source,(90,190,1220,910),zoom=lerp(0.86,1.05,ease(progress)))
    image.paste(plot,(90,190))
    draw.rounded_rectangle((1270,235,1830,815),radius=38,fill=(8,33,54,245),outline=ORANGE,width=7)
    draw.text((1550,325),"CURRENT MODEL",font=font(37,bold=True),fill=WHITE,anchor="mm")
    draw.text((1550,465),"~10",font=font(116,bold=True,rounded=True),fill=ORANGE,anchor="mm")
    draw.text((1550,565),"meV / atom",font=font(42,bold=True),fill=WHITE,anchor="mm")
    draw.text((1550,670),"energy-scale accuracy",font=font(29),fill=LIGHT_BLUE,anchor="mm")
    draw.rounded_rectangle((1305,720,1795,780),radius=16,fill=(245,128,37,230))
    draw.text((1550,750),"VALIDATE BEYOND PARITY",font=font(24,bold=True),fill=DARK_BROWN,anchor="mm")
    header(draw,"GRAPH-KERNEL ENERGY PREDICTION","FIRST RESULTS")
    footer(draw,"Phase stability + equation of state",progress)
    return image


def validation_scene(progress: float) -> Image.Image:
    image=Image.new("RGB",(WIDTH,HEIGHT),CREAM)
    draw=ImageDraw.Draw(image,"RGBA")
    source=Image.open(POSTER/"figures"/"bcc_magnetic_phonons.png").convert("RGB")
    plot=fit_image(source,(55,200,1285,915),zoom=0.92)
    image.paste(plot,(55,200))
    labels=[("FORCES",1400,300), ("EOS",1660,300), ("PHONONS",1400,520), ("FREE ENERGY",1660,520), ("MAGNETISM",1530,735)]
    for index,(label,x,y) in enumerate(labels):
        phase=(progress*5-index*0.13)%1
        color=ORANGE if phase<0.65 else LIGHT_BLUE
        draw.rounded_rectangle((x-120,y-58,x+120,y+58),radius=22,fill=(8,33,54,245),outline=color,width=6)
        draw.text((x,y),label,font=font(25,bold=True),fill=WHITE,anchor="mm")
    header(draw,"CHALLENGE THE MODEL WITH PHYSICS","VALIDATION BEYOND ENERGY")
    footer(draw,"Forces • EOS • phonons • free energies • magnetism",progress)
    return image


def poster_scene(progress: float) -> Image.Image:
    image=background((34,17,8),(92,42,14))
    draw=ImageDraw.Draw(image,"RGBA")
    source=Image.open(POSTER/"poster_preview.png").convert("RGB")
    zoom=lerp(0.80,0.88,ease(progress))
    target_w=int(1700*zoom/0.88); target_h=int(target_w*source.height/source.width)
    poster=source.resize((target_w,target_h),Image.Resampling.LANCZOS)
    x=(WIDTH-target_w)//2; y=165+(780-target_h)//2
    shadow=Image.new("RGBA",(target_w+50,target_h+50),(0,0,0,0))
    sd=ImageDraw.Draw(shadow); sd.rounded_rectangle((25,25,target_w+25,target_h+25),radius=18,fill=(0,0,0,150))
    shadow=shadow.filter(ImageFilter.GaussianBlur(16)); image.paste(shadow,(x-25,y-25),shadow)
    image.paste(poster,(x,y))
    draw.rounded_rectangle((370,735,1550,855),radius=32,fill=(8,33,54,245),outline=ORANGE,width=6)
    draw.text((960,790),"EXPLORE THE FULL POSTER",font=font(48,bold=True),fill=WHITE,anchor="mm")
    draw.text((960,845),"Supercomputing • Quantum simulations • Graph ML • Planetary science",font=font(25,bold=True),fill=ORANGE,anchor="mm")
    header(draw,"EARTH'S CORE THROUGH FIRST PRINCIPLES + MACHINE LEARNING","SC26 RESEARCH POSTERS")
    footer(draw,"Diego Juárez • The University of Texas at El Paso",progress)
    return image


SCENES = [
    (0.0, 10.0, earth_scene),
    (10.0, 20.0, hpc_scene),
    (20.0, 26.0, limits_scene),
    (26.0, 34.0, ml_scene),
    (34.0, 44.0, graph_scene),
    (44.0, 50.0, result_scene),
    (50.0, 56.0, validation_scene),
    (56.0, 60.0, poster_scene),
]


def render_at(time_s: float) -> Image.Image:
    for index,(start,end,renderer) in enumerate(SCENES):
        if time_s < end or index == len(SCENES)-1:
            progress=(time_s-start)/(end-start)
            current=renderer(progress)
            if index>0 and time_s-start<0.45:
                prev_start,prev_end,prev_renderer=SCENES[index-1]
                previous=prev_renderer(1.0)
                current=Image.blend(previous,current,ease((time_s-start)/0.45))
            return current
    return poster_scene(1.0)


def main():
    total_frames=FPS*DURATION
    writer=imageio.get_writer(
        OUTPUT,
        fps=FPS,
        codec="libx264",
        pixelformat="yuv420p",
        quality=8,
        macro_block_size=None,
        ffmpeg_log_level="warning",
        output_params=[
            "-vf", "scale=3840:2160:flags=lanczos,fps=30",
            "-movflags", "+faststart",
            "-metadata", "title=SC26 IronCoreMD Poster Blitz",
            "-metadata", "comment=Silent master for author-recorded narration",
        ],
    )
    try:
        for frame_index in range(total_frames):
            time_s=frame_index/FPS
            frame=np.asarray(render_at(time_s),dtype=np.uint8)
            writer.append_data(frame)
            if frame_index % (FPS*5)==0:
                print(f"rendered {frame_index/FPS:4.0f}/{DURATION} seconds",flush=True)
    finally:
        writer.close()
    print(OUTPUT)


if __name__ == "__main__":
    main()
