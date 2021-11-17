function [mask_arr] = test_seg(img_file_cell, rescale)
	% Arguments:
	%		img_file_cell -- one dimensional cell of image file path string. All the images have the same shape.
	% Returns:
	%		mask_arr -- array of masks for all the input images. The last dimension of the array equals to the length of the input image file cell.

	img_num = length(img_file_cell);

	ex_img_file = img_file_cell{1};
	ex_img = imread(ex_img_file);
	img_shape = size(ex_img);

	% scale_ratio = 0.5;
	scale_ratio = 1;
	% scale_ratio = 2;

	mask_arr = zeros(img_shape(1), img_shape(2), img_num);

	for i = 1:img_num
		img_file = img_file_cell{i};
		img = imread(img_file);

		if rescale
			img = imresize(img, scale_ratio);
		end

		%========================0: read image=====================================
		%==========================================================================
		sampleRGB=img;
		[height, width, channel]=size(sampleRGB);

		%========================1: color deconvolution============================
		%==========================================================================
		% SourceImage=sampleRGB;
		H_deinterlace = [0 1 0; 0 2 0; 0 1 0] ./4;

		sample_deinterlace = zeros(height, width, channel);
		for k=1:channel
			sample_deinterlace(:,:,k) = filter2(H_deinterlace,double(sampleRGB(:,:,k)),'same');
		end

		%=== Convert RGB intensity to optical density (absorbance)
		sampleRGB_OD = -log((sample_deinterlace+1)./256);

		% Construct color deconvolution matrix
		H2 = ones(10,10) ./ 100;
		sampleRGB_OD_Blur = zeros(height,width,channel);

		for k=1:channel
			sampleRGB_OD_Blur(:,:,k) = filter2(H2,sampleRGB_OD(:,:,k),'same');
		end

		% Standard values from literature
		He = [0.550 0.758 0.351]';
		Eo = [0.398 0.634 0.600]';
		Bg = [0.754 0.077 0.652]';

		% Create Deconvolution matrix
		M = [He/norm(He) Eo/norm(Eo) Bg/norm(Bg)];
		% D = inv(M);

		% Apply Color Deconvolution
		sampleHEB_OD = zeros(height, width, channel);
		for ki=1:height
			for ji=1:width
				RGB = reshape(sampleRGB_OD(ki,ji,:),channel,1);
				% HEB = D * RGB;
				HEB = M\RGB;

				sampleHEB_OD(ki,ji,1) = HEB(1);
				sampleHEB_OD(ki,ji,2) = HEB(2);
				sampleHEB_OD(ki,ji,3) = HEB(3);
			end
		end

		% Extract tumor cells that are stained by hematoxylin only
		hematoxylin = sampleHEB_OD(:,:,1);
		hematoxylin = (hematoxylin - 0.05) .* (hematoxylin > 0.05);
		hematoxylin = hematoxylin ./ max(max(hematoxylin));
		% h=fspecial('sobel');
		% hematoxylin_grad=sqrt(imfilter(hematoxylin,h,'replicate').^2+imfilter(hematoxylin,h','replicate').^2);


		%=================2: morphological operstions==============================
		%==========================================================================
		hImg=hematoxylin;

		se1=strel('disk',4);              %opening by reconstruction: remove small object
		hImgEro=imerode(hImg,se1);
		hImgEroRec=imreconstruct(hImgEro,hImg);

		se2=strel('disk',7);
		hImgDia=imdilate(hImgEroRec,se2);  %closing by reconstruction: remove small black spot
		hImgDiaRec=imreconstruct(imcomplement(hImgDia),imcomplement(hImgEroRec));
		hImgDiaRec=imcomplement(hImgDiaRec);
		hImgFill=imfill(hImgDiaRec,'holes');


		%=================3: Histogram Separation==================================
		%==========================================================================
		t = hImgFill;
		t1 = (t-min(min(t))) ./ (max(max(t)) - min(min(t))) .* 255;
		t1 = uint8(t1);
		[counts,binLocations] = imhist(t1);
		thresh_rm_small_mask = 50;
		% ind = find(counts>=thresh_rm_small_mask);
		ind = counts>=thresh_rm_small_mask;
		threshes = binLocations(ind);
		[h,w] = size(t1);
		masks = zeros(h,w,length(threshes));
		means = zeros(1,length(threshes));
		for ji = 1:length(threshes)
			mask = t1==threshes(ji);
			mask = bwareaopen(mask, thresh_rm_small_mask);
			mask = imfill(mask, 'holes');
			masks(:,:,ji) = mask;
			means(ji) = mean(mean(hematoxylin(mask)));
		end
		a0 = mean(mean(hematoxylin));
		target_mask_ind = find(means>=a0);
		final_mask = zeros(h,w);
		for ji = 1:length(target_mask_ind)
			final_mask = final_mask | masks(:,:,target_mask_ind(ji));
		end

		if rescale
			mask_arr(:,:,i) = imresize(center_select(final_mask), 1/scale_ratio);
		else
			mask_arr(:,:,i) = center_select(final_mask);
		end
	end


	% se3=strel('disk',1);
	% final_mask_1 = imclose(final_mask, se3);
	% figure; imshow(final_mask_1); title('final mask 1');
	% final_mask_2 = imopen(final_mask_1, se3);
	% figure; imshow(final_mask_2); title('final mask 2');

	% %****************** Watershed implementation from MATLAB Documentation ******
	% bw = final_mask_1;
	% D = bwdist(~bw);
	% D = -D;
	% L = watershed(D);
	% %***************** This is the segmented objects' labels ***********
	% L(~bw) = 0;

	% %***************** Display the resulting label matrix as an RGB image *******
	% %***************** Comment when you do not need it *****************
	% rgb = label2rgb(L,'jet',[.5 .5 .5]);
	% figure; imshow(rgb)
end











