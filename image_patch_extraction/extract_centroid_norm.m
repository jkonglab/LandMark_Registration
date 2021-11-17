function [] = extract_centroid_norm(img_file, magnification, save_name)

	% Input:
	% img_file			- string, image file to be processed.
	% magnification		- integer.
	% save_name			- string, file to store centroid mask image.

	img = imread(img_file);
	img_norm = normalizeStaining(img);
	mask = segNuclei(img_norm, magnification);
	imwrite(mask, save_name);

end










