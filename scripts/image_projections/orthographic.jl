println("Orthographic")

dest_proj = Orthographic(;long0=0.0, lat0=30.0)
back_long0 = dest_proj.long0 <= 0 ? dest_proj.long0 + 180 : dest_proj.long0 - 180
back_lat0 = -dest_proj.lat0
dest_proj_back = Orthographic(;long0=back_long0, lat0=back_lat0)

width_dest, height_dest = 500, 500
strip_width = 25

height, width = size(img)
r = dest_proj.radius
dest_affine = affine_from_bounds(-r, -r, r, r, width_dest, height_dest);

out_img_front = reproject_warp(
    img, (height_dest, width_dest), src_proj, dest_proj, src_affine, dest_affine; method=bilinear)
out_img_back = reproject_warp(
    img, (height_dest, width_dest), src_proj, dest_proj_back, src_affine, dest_affine; method=bilinear)
strip = zeros(RGB, height_dest, strip_width)
out_img = hcat(out_img_front, strip, out_img_back)
out_img = map(clamp01nan, out_img)

figure_size = (max_figure_size, ceil(Int, max_figure_size/2))

canvas = plot(
    out_img, 
    xlims=(1, 2 * width_dest),
    ylims=(1, height_dest),
    size=figure_size, 
    ticks=:none,
    )

gridlines = make_grid_lines(
    dest_affine, 
    x -> project(dest_proj, x, clip=true),
    -180:30.0:180.0, 
    -90.0:30.0:90;
    num_points=1000
);
for line in gridlines
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs, ys, color=:white, label="", alpha=alpha)
end
gridlines = make_grid_lines(
    dest_affine, x -> project(dest_proj_back, x, clip=true), 
    -180:30.0:180.0,
    -90.0:30.0:90;
    num_points=1000
    );
for line in gridlines
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs .+ strip_width .+ width_dest, ys, color=:white, label="", alpha=alpha)
end

output_path = joinpath(output_dir, "orthographic.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")