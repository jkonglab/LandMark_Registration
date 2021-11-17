function [out_mask] = center_select(in_mask)

	% Select the target near the center of the input mask. Then fill the target's holes and try to make it to be more convex. Since most nucleus are of the convex shape.
	%
	% Argument:
	%		in_mask -- input binary mask image.
	%
	% Returns:
	%		out_mask -- output binary mask image.

	mask_size = size(in_mask);
	H = mask_size(1);
	W = mask_size(2);
	center_mask = zeros(H, W);
	out_mask = zeros(H, W);
	center_mask(round(H/2)-5:round(H/2)+5, round(W/2)-5:round(W/2)+5) = 1;

	label = bwlabel(in_mask);

	intersect_label = unique(label(logical(center_mask)));
	intersect_label = intersect_label(intersect_label~=0);

	if ~ isempty(intersect_label)
		for li = 1:length(intersect_label)
			out_mask = out_mask | label==intersect_label(li);
		end

		% out_mask = imfill(out_mask, 8, 'holes');

		% se = strel('disk',5);
		% se = strel('disk',10);
		se = strel('disk',20);
		out_mask = imclose(out_mask, se);

		out_mask = imfill(out_mask, 8, 'holes');
	end

end

