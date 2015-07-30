function [] = convertDataTo32Tiff( filename, data )

result = Tiff(filename,'w');
data = im2single(data);
tagstruct.ImageLength = size(data,1);
tagstruct.ImageWidth = size(data,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
tagstruct.SampleFormat = 3;
result.setTag(tagstruct);
result.write(data);
result.close();

end 