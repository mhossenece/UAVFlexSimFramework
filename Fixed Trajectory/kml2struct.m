% Edited by Md Sharif Hossen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function kmlStruct = kml2struct(kmlFile)
    % Import a .kml file as a vector array of shapefile structs with fields:
    % Geometry, Name, Description, Lon, Lat, and BoundingBox.
    % The function supports points, lines, and polygons.

    % Open the file
    [FID, msg] = fopen(kmlFile, 'rt');
    if FID < 0
        error(msg)
    end
    txt = fread(FID, 'uint8=>char')';
    fclose(FID);

    % Find and parse the schema
    expr = '<Schema .+?>[\S\s]*?</Schema>';
    schemaString = regexp(txt, expr, 'match');
    if ~isempty(schemaString)
        expr = '<SimpleField .+?>[\S\s]*?</SimpleField>';
        schemaFields = regexp(schemaString{1}, expr, 'match');
        schemaFields = regexprep(schemaFields, '<SimpleField name="', '');
        schemaFields = regexprep(schemaFields, '".*', '');
    else
        schemaFields = {};
    end

    % Find the placemarks and put them in an array
    expr = '<Placemark.+?>.+?</Placemark>';
    objectStrings = regexp(txt, expr, 'match');
    Nos = length(objectStrings);

    % Initialize output struct array
    kmlStruct = struct('Geometry', '', 'Name', '', 'Description', '', 'Lon', [], 'Lat', [], 'BoundingBox', []);

    % Loop through each placemark and extract data
    for ii = 1:Nos
        % Find Object Name Field
        bucket = regexp(objectStrings{ii}, '<name.*?>.*?</name>', 'match');
        if isempty(bucket)
            name = 'undefined';
        else
            name = regexprep(bucket{1}, '<name.*?>\s*', '');
            name = regexprep(name, '\s*</name>', '');
        end
        kmlStruct(ii).Name = name;

        % Find Object Description Field
        bucket = regexp(objectStrings{ii}, '<description.*?>.+?</description>', 'match');
        if isempty(bucket)
            desc = '';
        else
            desc = regexprep(bucket{1}, '<description.*?>\s*', '');
            desc = regexprep(desc, '\s*</description>', '');
        end
        kmlStruct(ii).Description = desc;

        geom = 0;
        % Identify Object Type
        if ~isempty(regexp(objectStrings{ii}, '<Point', 'once'))
            geom = 1;
        elseif ~isempty(regexp(objectStrings{ii}, '<LineString', 'once'))
            geom = 2;
        elseif ~isempty(regexp(objectStrings{ii}, '<Polygon', 'once'))
            geom = 3;
        end

        switch geom
            case 1
                geometry = 'Point';
            case 2
                geometry = 'Line';
            case 3
                geometry = 'Polygon';
            otherwise
                geometry = '';
        end
        kmlStruct(ii).Geometry = geometry;

        % Find Coordinate Field
        bucket = regexp(objectStrings{ii}, '<coordinates.*?>.*?</coordinates>', 'match');
        if ~isempty(bucket)
            coordStr = regexprep(bucket{1}, '<coordinates.*?>(\s+)*', '');
            coordStr = regexprep(coordStr, '(\s+)*</coordinates>', '');
            coordMat = str2double(regexp(coordStr, '[,\s]+', 'split'));
            coordMat = reshape(coordMat, 3, [])';

            Lon = coordMat(:, 1);
            Lat = coordMat(:, 2);

            if geom == 3
                Lon = [Lon; NaN];
                Lat = [Lat; NaN];
            end

            kmlStruct(ii).Lon = Lon;
            kmlStruct(ii).Lat = Lat;
            kmlStruct(ii).BoundingBox = [min(Lon) min(Lat); max(Lon) max(Lat)];
        else
            kmlStruct(ii).Lon = [];
            kmlStruct(ii).Lat = [];
            kmlStruct(ii).BoundingBox = [];
        end

        % Schema Object Fields
        for jj = 1:length(schemaFields)
            expr = strcat('<SimpleData.*?name="', schemaFields{jj}, '".*?>.*?</SimpleData>');
            bucket = regexp(objectStrings{ii}, expr, 'match');
            if isempty(bucket)
                name = 'undefined';
            else
                expr = strcat('<SimpleData .*?name="', schemaFields{jj}, '".*?>\s*');
                name = regexprep(bucket{1}, expr, '');
                name = regexprep(name, '\s*</SimpleData>', '');
            end
            kmlStruct(ii).(schemaFields{jj}) = name;
        end
    end
end
