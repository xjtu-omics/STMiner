def cut_image(image):
    height, width = image.shape[:2]

    top = 0
    bottom = height - 1
    left = 0
    right = width - 1

    while top < height:
        if image[top].any():
            break
        top += 1

    while bottom >= 0:
        if image[bottom].any():
            break
        bottom -= 1

    while left < width:
        if image[:, left].any():
            break
        left += 1

    while right >= 0:
        if image[:, right].any():
            break
        right -= 1

    cropped_image = image[top:bottom + 1, left:right + 1]
    return cropped_image
