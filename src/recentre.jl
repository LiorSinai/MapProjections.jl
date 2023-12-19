"""
    recentre(image, affine, long0)

Recentre an image at `long0`.
"""
function recentre(img::AbstractMatrix, affine::AffineTransform, long0::AbstractFloat)
    height, width = size(img)
    j0, i0 = inv(affine)((long0, 0.0))
    left = floor(Int, j0 - width / 2)
    right = ceil(Int, j0 + width / 2)
    if right >= width
        js = vcat(left:width, 1:(left-1))
    elseif left < width
        js = vcat(right:width, 1:(right-1))
    end
    img[:, js]
end
