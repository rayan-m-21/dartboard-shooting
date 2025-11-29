import cv2
import params as params
import numpy as np

# Conversion scale: mm to pixels
IMG_SIZE = 2000
CENTER = (IMG_SIZE // 2, IMG_SIZE // 2)
MM_PER_PX = 2 * (IMG_SIZE / 1000)

# Convert to pixels
radii_px = {k: int(v * MM_PER_PX) for k, v in params.radii_mm.items()}

# Colors (BGR format)
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (0, 0, 255)
GREEN = (0, 255, 0)
DARK_GREEN = (0, 100, 0)
DARK_RED = (0, 0, 139)

def draw_numbers_on_board(img):
    """Draws numbers 1–20 around the dartboard."""
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_scale = int(10 * (IMG_SIZE / 10000))
    thickness = int(30 * (IMG_SIZE / 10000))

    number_radius_mm = params.radii_mm['double_outer'] + 15  # place numbers just outside the edge
    for i, number in enumerate(params.start_angles):
        angle_deg = params.start_angles[number] + 9  # add half of segment (18° / 2)
        angle_rad = np.deg2rad(angle_deg)

        # Compute position
        x_mm = number_radius_mm * np.cos(angle_rad)
        y_mm = number_radius_mm * np.sin(angle_rad)
        x_px, y_px = mm_to_pix((x_mm, y_mm))

        text = str(number)
        text_size = cv2.getTextSize(text, font, font_scale, thickness)[0]
        text_x = x_px - text_size[0] // 2
        text_y = y_px + text_size[1] // 2

        cv2.putText(img, text, (text_x, text_y), font, font_scale, BLACK, thickness, cv2.LINE_AA)

def get_dartboard_image():
    """
    Returns the openCV2 image for the blank dartboard with the
    above params
    """
    # Create blank image
    img = np.ones((IMG_SIZE, IMG_SIZE, 3), dtype=np.uint8) * 255

    # Draw outer bull (green) and inner bull (red)
    cv2.circle(img, CENTER, radii_px['outer_bull'], GREEN, -1)
    cv2.circle(img, CENTER, radii_px['inner_bull'], RED, -1)

    # Draw 20 segments and color them
    for i in range(20):
        start_angle = -9 + i * 18
        end_angle = start_angle + 18

        # Alternate black and white for base segments
        base_color = WHITE if i % 2 == 0 else BLACK
        cv2.ellipse(img, CENTER, (radii_px['double_inner'], radii_px['double_inner']),
                    0, start_angle, end_angle, base_color, -1)

        # Double ring: alternate red/green
        ring_color = DARK_GREEN if i % 2 == 0 else DARK_RED
        cv2.ellipse(img, CENTER, (radii_px['double_outer'], radii_px['double_outer']),
                    0, start_angle, end_angle, ring_color, -1)
        cv2.ellipse(img, CENTER, (radii_px['double_inner'], radii_px['double_inner']),
                    0, start_angle, end_angle, base_color, -1)

        # Triple ring: alternate red/green
        ring_color = DARK_GREEN if i % 2 == 0 else DARK_RED
        cv2.ellipse(img, CENTER, (radii_px['triple_outer'], radii_px['triple_outer']),
                    0, start_angle, end_angle, ring_color, -1)
        cv2.ellipse(img, CENTER, (radii_px['triple_inner'], radii_px['triple_inner']),
                    0, start_angle, end_angle, base_color, -1)

    # Draw inner circle
    cv2.circle(img, CENTER, radii_px['outer_bull'], DARK_GREEN, -1)
    cv2.circle(img, CENTER, radii_px['inner_bull'], RED, -1)

    # Draw circular outlines for clarity
    for key in ['board_edge', 'double_outer', 'double_inner', 'triple_outer', 'triple_inner', 'outer_bull', 'inner_bull']:
        cv2.circle(img, CENTER, radii_px[key], BLACK, int(2 * (IMG_SIZE / 10000)))

    # draw letters
    draw_numbers_on_board(img)

    # return image
    return img

def draw_sequence(img, points):
    """
    Draws a list of points with a red-to-blue color gradient on the image.
    """
    num_points = len(points)
    for i, (x, y) in enumerate(points):
        # Interpolate between red (0,0,255) and blue (255,0,0)
        alpha = i / max(num_points - 1, 1)
        b = int(255 * alpha)
        r = int(255 * (1 - alpha))
        color = (b, 0, r)  # BGR format
        cv2.circle(img, (int(x), int(y)), 15, color, -1)

    return img

def mm_to_pix(cartesian_coord):
    """
    Converts (x_mm, y_mm) in mm (Cartesian, origin at center)
    to (x_px, y_px) in image coordinates (origin top-left).
    """
    x_mm, y_mm = cartesian_coord
    x_px = int(CENTER[0] + x_mm * MM_PER_PX)
    y_px = int(CENTER[1] - y_mm * MM_PER_PX)  # y is flipped
    return (x_px, y_px)