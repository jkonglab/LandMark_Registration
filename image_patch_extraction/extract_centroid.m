function mask = extract_centroid(img, magnification)

	% Input:
	% img_file			- string, image file to be processed.
	% magnification		- integer.
	% save_name			- string, file to store centroid mask image.

	%img = imread(img_file);
	mask = segNuclei(img, magnification);
	%imwrite(mask, save_name);

end










