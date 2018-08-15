# LineThruPoly2
whew
function [ dist_total, dist_straight_line,weighted_straight_line] ...
    = LineThruPoly2(polygons, begindex, endex, node_from_index, node_to_index, node_from_location, node_to_location, polycost,resolution, ...
    boundary,tolerance, lakejump, proximity_to_boundary,shortestpath_bool, shortestpath_lb, shortestpath_ub, reference_latitude)
tic
%function will accept vector of transmission lines from, and to, determine
%their Euclidean distance in the plane using weight functions to describe
%passing through "high cost" areas such as lakes or mountains and penalize
%accordingly The function will determine if it
%is cheaper to go 'around' the polygon or through it using a shortest path
%algorithm when it is determined to be necessary via linear correlation.
%the function inputs are node indices corresponding to connection, node
%location in LONGITUDE then LATITUDE (CORRESPONDING TO X AND Y), polycost which is a linear coefficient
%describing the 'cost' of going through a polygon, and resolution creates a
%resolution-by-resolution mesh when the shortest path algorithm is called
%endex stores the ending indices for each of the shape vectors to save time
%and effort.

%   output is a vector of total distance of each transmission line determined by  
%   euclidean distance and weight of the blocking polygons
%   INPUT LATITUDE THEN LONGITUDE FOR 
%   NODE LOCATION. polygonx is longitude, polygony is latitude. (yes i know
%   that is confusing, that's how the NUTS data was given and how i
%   received location of node data).
dist_total = zeros(length(node_from_index),1); %dist_total returns the total distance that a connection should have to travel, minimized using shortest path
dist_straight_line = zeros(length(node_from_index),1); %dist_straight_line returns simply the straight-line distance between connections without considering blocking polygons
weighted_straight_line = zeros(length(node_from_index),1); %weighted_straihgt_line returns the straight-line distance between connections BUT it considers any and all introduced blocking and resident polygons as punishments. 
if reference_latitude == 6371
    R = reference_latitude;
else
    reference_latitude = reference_latitude*pi/180;
    R = 6378 - 21*sin(reference_latitude);
end
node_from = [node_from_index, node_from_location];
node_to = [node_to_index, node_to_location];
new_node_from = zeros(length(node_from_index),3);
for i = 1:length(node_from_index) %switch all node_froms to be the minimum node index in the pair
    if node_from(i,1) > node_to(i,1)
        new_node_from(i,:) = node_to(i,:);
        node_to(i,:) = node_from(i,:);
        node_from(i,:) = new_node_from(i,:);
    end
end
connections = [node_from, node_to];
%connections = sortrows(connections,1); %will sort rows based on node_froms        
all_nodes = unique([node_from; node_to],'rows'); %will return the unique rows and reject the ones that are identical i.e. same index and location
all_nodes_x = all_nodes(:,2); %collects longitudes for y data
all_nodes_y = all_nodes(:,3); %collects latitudes for x data
row_polygon = length(begindex);
index_bus_in_poly = zeros(length(all_nodes),row_polygon);

[num_connections,~] = size(node_from_location);
[row_nodes,~] = size(all_nodes);
%firstly, can find nodes in polygon
for i = 1:row_nodes
    for j = 1:row_polygon
        if inpolygon(all_nodes_x(i),all_nodes_y(i),polygons(begindex(j):endex(j),1), polygons(begindex(j):endex(j),2))
            index_bus_in_poly(i,j) = 1;
        end
    end
end


    
still_in_poly = ones(size(index_bus_in_poly));
index_bus_in_poly = sparse(index_bus_in_poly); %returns the node index i in polygon j


[I,J] = find(index_bus_in_poly); %finding the I nodes and their corresponding J polygons. ith row corresponds to the ith node, jth column corresp to jth polygon

%need to determine which buses in poly are not actually in the poly
%because of pre-defined tolerance 

for i  = 1:length(I) %over all nodes determined to be in polygons 
    if tolerance(J(i)) ~= 0 %if there is a nonzero tolerance for the polygon corresponding to this node
        for j = begindex(J(i)):endex(J(i)) % over the length of each of the vectors of lat and long corresponding to each polygon
            if ~isnan(polygons(j,2)) %if it's not a nan we can check to see what the distance is between the points
                dist_from_boundary = deg2km(sqrt((all_nodes_x(I(i)) - polygons(j,1))^2 + ...
                (all_nodes_y(I(i)) - polygons(j,2))^2),R);
                if dist_from_boundary < tolerance(J(i)) %if it's determined that the distance is within the tolerance, then the node will no longer be considered as in the poly
                    still_in_poly(I(i),J(i)) = 0; %still in poly is an identically sized matrix 
                    break %break the while loop as we have already established that we are within the tolerance
                end
            end
            
        end    
    end
    
end
        

%selecting connections to analyze. Find the square that best defines each polygon and
%the line that defines the connection. If any of the four corners of the
%square are independent from the others with respect to which side of the line they are on
%, then there must be an intersection. 
check_matrix = zeros(num_connections,row_polygon);
for i = 1:num_connections
    for j = 1:row_polygon
        if i == 343 && j == 147
            doot = 1;
        end
        if boundary(j) == 0 %if it's a bad polygon we can check the square thing
            max_bound_x = max(polygons(begindex(j):endex(j),1));
            max_bound_y = max(polygons(begindex(j):endex(j),2));
            min_bound_x = min(polygons(begindex(j):endex(j),1));
            min_bound_y = min(polygons(begindex(j):endex(j),2));
            top_right = [max_bound_x, max_bound_y];
            top_left = [min_bound_x, max_bound_y];
            bottom_left = [min_bound_x, min_bound_y];
            bottom_right = [max_bound_x, min_bound_y];
            four_corners = [top_right; top_left; bottom_left; bottom_right];
            line_connector = [linspace(connections(i,2),connections(i,5)); linspace(connections(i,3)...
                ,connections(i,6))];
            miner1 = 1e+6;
            miner2 = 1e+6;
            for k = 1:size(line_connector,2) %find which point in the vectorization that the min and max lie at
                checkmax = abs(line_connector(1,k)-max_bound_x);
                checkmin = abs(line_connector(1,k)-min_bound_x);
                if checkmax < miner1
                    miner1 = checkmax;
                    maxdex = k;
                end
                if checkmin < miner2
                    miner2 = checkmin;
                    mindex = k;
                end
                
            end
            above_or_below = zeros(size(four_corners,1),1);
            for k = 1:size(four_corners,1)
                if four_corners(k,1) == max_bound_x
                    if line_connector(2,maxdex) > four_corners(k,2) %if the line is above the corner
                        above_or_below(k) = 0; %0 means below
                    else
                        above_or_below(k) = 1; %1 means above 
                    end
                else %otherwise it is a minimum, perform same check
                    if line_connector(2,mindex) > four_corners(k,2)
                        above_or_below(k) = 0;
                    else
                        above_or_below(k) = 1;
                    end
                end
            end
            if length(unique(above_or_below)) == 2 %if the vector has both 1s and 0s
                check_matrix(i,j) = 1; %submit the connection to the check matrix. 
            end
        else %it's a good polygon then we check the proximity to the boundary to determine whether we want to submit it to the check matrix
            
            minimum = proximity_to_boundary;
            for k = begindex(j):endex(j) %over all the points corresp to the good polygon
                dist_from_poly = deg2km(sqrt((connections(i,2)-polygons(k,1))^2 + (connections(i,3)-polygons(k,2))^2),R);
                if dist_from_poly < minimum
                    minimum = dist_from_poly;
                    if minimum < proximity_to_boundary
                        check_matrix(i,j) = 1;
                        break %once we've established we want to check the ith connection against the jth polygon, break the k loop
                    end
                end
            end
        end
    end
end
check_matrix = sparse(check_matrix);
[M,N] = find(check_matrix); %M is the set of connections to be checked against the N polygons with which they could possibly be incident (m,n) is one such pair

blocking_polys = zeros(size(check_matrix)); %confirms that a poly n will block a connection m
x_intersect = zeros([size(check_matrix) 200]); %x_intersect will store in a 3x3 matrix connection m's crossings with polygon n
y_intersect = zeros([size(check_matrix) 200]); %actually i
%feel like i don't need this because you're just gonna evaluate the
%shortest path anyway, although maybe it would save memory
for i = 1:length(M)
    [x_cross, y_cross] = polyxpoly(linspace(connections(M(i),2), connections(M(i),5),10),linspace(connections(M(i),3), connections(M(i),6),10), ...
        polygons(begindex(N(i)):endex(N(i)),1), polygons(begindex(N(i)):endex(N(i)),2));
    if ~isempty(x_cross) %if there are crossings, then we indicate that in the blocking_polys matrix 
        blocking_polys(M(i),N(i)) = 1; %confirms that checked connection M(i) is intersected by a corresponding polygon N(i)
        x_intersect(M(i),N(i),1:length(x_cross)) = x_cross;
        y_intersect(M(i),N(i),1:length(y_cross)) = y_cross;
    end 
end
blocking_polys = sparse(blocking_polys);
[P,Q] = find(blocking_polys); %stores the indices of a connection p blocked by a poly q 
i = 1;
while 1 %chain of if statements to evaluate the properites of the connection with respect to location. i corresponds to the ith connection
    if i > num_connections
        break
    end
    lat2 = connections(i,6)*pi/180;
    lat1 = connections(i,3)*pi/180;
    long2 = connections(i,5)*pi/180;
    long1 = connections(i,2)*pi/180;
    diff_long = long2-long1; %finding the difference in longitude. 
    diff_lat = lat2-lat1; %finding the difference in latitude.
    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2;
    c = 2*asin(min(1,sqrt(a)));
    dist_straight_line(i) = c * R;
    connection_blocked = find(P==i); %finds whether connection i is an element of the set of blocked connections P, and also returning the id of the blocking polygon
    
    %firstly we can collect the values of the blocking poly intersections
    %with find
    [~,blocked_indexx,x_blockers] = find(x_intersect(i,Q(connection_blocked),:));
    [~,~,y_blockers] = find(y_intersect(i,Q(connection_blocked),:));
    blocked_indexx = mod(blocked_indexx, length(connection_blocked)); %split the values by their associated blocking polygon
    blocked_indexx(~blocked_indexx) = length(connection_blocked); %replace zero mods with the last number. Now you have the indices in connection_blocked corresponding to the associated polygon.
    %next we can determine if the principal or terminal nodes are in
    %polygons
    if ismember(connections(i,1),I) %if the first node is contained within a polygon
        j = J(I==connections(i,1)); %returning the value of all polygons that contain the principal node
        j = unique(j);
    end
    if ismember(connections(i,4),I) %if the terminal node is also in a polygon
        k = J(I==connections(i,4)); %store the values of all polygons that contain the principal node
        k = unique(k);
    end
    %if there are no blocking polygons then maybe we can ignore them
    %entirely and just do a straight line calculation. 
    if ~isempty(x_blockers) %if x_blockers is not empty, then there are blocking polys (which will always be the case for a scenario where a node is inside a polygon unless the terminal node is in the same polygon but wait that's not even true necessarily. 
        sequencex = [x_blockers; connections(i,5); connections(i,2)];
        sequencey = [y_blockers; connections(i,6); connections(i,3)];
        start_index = length(sequencex);
                
        sequence = sequencer(sequencex, sequencey, start_index); %sequencer returns the indices of the way the line is connected, sorted by blockages on the path from principal to terminal.
                %calculate weighted straight line distance by piecing
                %together along the sequence. Evaluate whether or not to
                %punish based on lakejump. Also need to evaluate whether to
                %punish the nodes inside polygons based on tolerance.
                %Firstly collect all the values of the middle of each of
                %the connections, and then check to see how many of them
                %are inside the polygons that also contain the node. This
                %will help us decide how to punish the line based on
                %tolerance. 
        long1 = sequencex(sequence(1:end-1));
        long2 = sequencex(sequence(2:end));
        lat1 = sequencey(sequence(1:end-1));
        lat2 = sequencey(sequence(2:end));
        mean_long = (long1 + long2)./2; %find the middle of each line segment in the connection
        mean_lat = (lat1 + lat2)./2;
        connection_contained = zeros(length(mean_long), length(connection_blocked));
        for count = 1:length(long2) %now need to check which polygons contain each of the means. 
              for count2 = 1:length(connection_blocked)
                   if inpolygon(long2(count),lat2(count),polygons(begindex(Q(connection_blocked(count2))):endex(Q(connection_blocked(count2))),1),...
                                polygons(begindex(Q(connection_blocked(count2))):endex(Q(connection_blocked(count2))),2)) && ...
                                inpolygon(mean_long(count),mean_lat(count),polygons(begindex(Q(connection_blocked(count2))):endex(Q(connection_blocked(count2))),1),...
                                polygons(begindex(Q(connection_blocked(count2))):endex(Q(connection_blocked(count2))),2))%if the end of the line and the middle of the line are both in or on the shape, then we can logi
                            
                            connection_contained(count,count2) = 1;
                   end
              end
        end
                %shared_home = connection_contained(:,j); %return all the columns corresponding to the polygons that contain the principal node
                %shared_away = connection_contained(:,k); %return all the columns corresponding to the polygons that contain the terminal node
                %[home_connection, home_poly] = find(shared_home); %get the values of all the connections that are inside the polygon(s) that contain the princip node
                %[away_connection, away_poly] = find(shared_away);
                %determine which of the principal node's polygons to
                %punish. 
        num_segments = length(mean_long);
                
        for count = 1:num_segments %here we will add up the segments and their associated penalties, considering variables like tolerance and lakejump
            polys = find(connection_contained(count,:)); %returns the index of all the polys that contain this segment.
            if any(ismember(Q(connection_blocked(polys)), j)) %if any of the polys that contain this segment are the same as the principal node, then we want to perform tolerance check. 
                which_ones = ismember(Q(connection_blocked(polys)),j);
                for count2 = 1:length(which_ones)
                    if which_ones(count2) %if this one is a member of the polygon, see if the far part is farther away than the tolerance
                        distance = deg2km(sqrt((connections(i,2)-long2(count))^2 + (connections(i,3)-long2(count))^2),R);
                         %used for determining the exit point of the line from the principal polygon into free space which is good for starting the mesh. Although do you really wanna....
                        diff_long = long2(count) - long1(count);
                        diff_lat = lat2(count) - lat1(count);
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2;
                        c = 2*asin(min(1,sqrt(a)));
                        d = c * R;
                        if distance > tolerance(Q(connection_blocked(polys(count)))) %if the distance is less than the tolerance, we don't have to punish this segment 
                            weighted_straight_line(i) = weighted_straight_line(i) + d*(polycost(Q(connection_blocked(polys(count))))-1);
                        end
                        inside_index_principal = count2; %storing the last index associated with the line's exit from the polygons that bound the principal node. 
                    end
                end
            elseif any(ismember(Q(connection_blocked(polys)),k)) %it is a member of any of the polygons that contain the terminal node
                which_ones = ismember(Q(connection_blocked(polys)),k); %logical vector stating which of the polygons that this line is a part of also contain the terminal node
                for count2 = length(which_ones):1
                    if which_ones(count2) %if this one is a member of a polygon containing the terminal node
                        distance = deg2km(sqrt((connections(i,5)-long1(count))^2 + (connections(i,6)-long1(count))^2),R);
                        
                        diff_long = long2(count) - long1(count);
                        diff_lat = lat2(count) - lat1(count);
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2;
                        c = 2*asin(min(1,sqrt(a)));
                        d = c*R;
                        if distance > tolerance(Q(connection_blocked(polys(count))))
                            weighted_straight_line(i) = weighted_straight_line(i) + d*(polycost(Q(connection_blocked(polys(count))))-1);
                        end
                        inside_index_terminal = count2;
                                
                    end
                end
            else %it is not a member of any of the polygons that contain the principal or terminal node. Do lakejump. 
                diff_long = long2(count) - long1(count);
                diff_lat = lat2(count) - lat1(count);
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2;
                c = 2*asin(min(1,sqrt(a)));
                d = c * R;
                if ~lakejump || d > resolution %if we're not lakejumping or the distance is greater than the resolution, then we need to dub out and punish
                    weighted_straight_line(i) = weighted_straight_line(i) + d*(polycost(Q(connection_blocked(polys(count))))-1);
                    %we can lakejump. not possible to remove the penalty so
                    %we have to decide something else... maybe make a
                    %square around it with an inverse penalty so they
                    %cancel and it keeps on through :) 
                    
                    square_max_long = max([long1(count) long2(count)]);
                    square_min_long = min([long1(count) long2(count)]);
                    square_max_lat = max([lat1(count) lat2(count)]);
                    square_min_lat = min([lat1(count) lat2(count)]);
                    squaremaxlong = square_max_long + resolution;
                    squareminlong = square_min_long - resolution;
                    squaremaxlat = square_max_lat + resolution;
                    squareminlat = square_min_lat - resolution;
                    negate_square_x = [squaremaxlong; squaremaxlong; squareminlong; squareminlong; squaremaxlong]; %a bounding square for the cancellation of the penalty for passing through this strip 
                    negate_square_y = [squareminlat; squaremaxlat; squaremaxlat; squareminlat; squareminlat];
                else 
                end
            end
            weighted_straight_line(i) = weighted_straight_line(i) + dist_straight_line(i); %already added weighted parts of segments, this is the rest of that bitch. 
        end
            %now need to implement shortest path algorithm call. Need to
            %consider blockages. Do you want point of exit of the polygons
            %or are you just gonna throw everything in . I'm thinking the
            %latter option to avoid the headache.
    else %there are no blocking polys. just add em up.
        weighted_straight_line(i) = dist_straight_line(i);
        dist_total(i) = dist_straight_line(i);
    end
                %for each column, must punish the maximum value in each
                %column according to the polycost because those are the
                %ones that stretch the farthest. But also need to ignore
                %the summation of everything else. For example if you want
                %to ignore the first polygon but not the second, only want
                %the summation of the line segments that contain the second
                %punishment which should be the first two line segments. Be
                %careful not to double count the summations. 
                
                
                
                
                
                
                
                for count = 1:length(sequence)-1
                    if count == 1 %we are at the principal node, inside the first polygon. We can determine whether or not to punish based on still_in_poly
                        long1 = sequencex(sequence(count));
                        long2 = sequencex(sequence(count+1));
                        lat1 = sequencey(sequence(count));
                        lat2 = sequencey(sequence(count+1));
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2;
                        c = 2*asin(min(1,sqrt(a)));
                        d = c*R;
                        if still_in_poly(connections(i,1),j) %if the connection is still going to be punished
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(j);
                        else %it's not still in the poly and we don't punish
                            weighted_straight_line(i) = weighted_straight_line(i) + d;
                        end
                    elseif count == length(sequence)-1 %we are at the final counter thing and need to evaluate as above
                        long1 = sequencex(sequence(count));
                        long2 = sequencex(sequence(count+1));
                        lat1 = sequencey(sequence(count));
                        lat2 = sequencey(sequence(count+1));
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2;
                        c = 2*asin(min(1,sqrt(a)));
                        d = c*R;
                        if still_in_poly(connections(i,4),k) %if the terminal node is still to be punished
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(k);
                        else %don't punish the terminal node. 
                            weighted_straight_line(i) = weighted_straight_line(i) + d;
                        end
                    else % we are at intermediates and need to perform ~analysis~ to determine whether they are crossing anything. 
                        long1 = sequencex(sequence(count));
                        long2 = sequencex(sequence(count+1));
                        lat1 = sequencey(sequence(count));
                        lat2 = sequencey(sequence(count+1));
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2;
                        c = 2*asin(min(1,sqrt(a)));
                        d = c*R;
                        long_midpoint = mean([long1 long2]); %if the midpoint is inside or outside a polygon you'll need to weight them separately. Also need to consider embedded polygons. 
                        lat_midpoint = mean([lat1 lat2]);
                        %need the polygon associated with these midpoints.
                        %Access via
                        %Q(connection_blocked(blocked_indexx(sequence(count))))
                        polycheck = Q(connection_blocked(blocked_indexx(sequence))); %check if midpoint is inside first polygon
                        punish = zeros(size(connection_blocked)); %zeros. Will add polycost if it's decided the line must be punished. Gonna check all of the connection blockers because i'm not smart enough to do it any other way. 
                        for count2 = 1:length(polycheck)
                            if inpolygon(long_midpoint, lat_midpoint, polygons(begindex(polycheck(count2)):endex(polycheck(count2)),1), ...
                                polygons(begindex(polycheck(count2)):endex(polycheck(count2)),2)) %if inside the first polygon, punish the connection
                            
                            %punish the line according to the polycost
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(polycheck(count2));
                            else %it's not in that particular polygon, don't penalize it
                            end
                        end
                        
                    end
                end
            else %else there are no blocking polygon intersections which means the nodes are in the same polygon
                weighted_straight_line(i) = dist_straight_line*polycost(j);
                dist_total(i) = weighted_straight_line(i);
                
    end
            
end
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if ismember(connections(i,1),I) %if the first node is contained within a polygon
        j = J(find(I==connections(i,1),1)); %corresp polygon to node´s location. only need to find first instance because the corresp poly will always be the same for the same principal node
            if ismember(connections(i,4),I) %if the terminal node is in a polygon
                k = J(find(I==connections(i,4),1)); %find the corresp polygon for the terminal node
                if j == k
                    %if the nodes are in the same poly, then distance is simply
                    %between them times the polycost SCENARIO 1
                   
                            dist_total(i) = dist_straight_line(i)*polycost(j); %total distance of connection
                            weighted_straight_line(i) = dist_total(i); %the weighted straight line distance is the total distance of the connection. 
                        
                
                
                
                elseif length(connection_blocked)>2 %if there exists at least one blocking poly (not including the poly containing princ. and term. nodes) corresponding to this connection
                    %SCENARIO 2
                    poly_in = Q==j; %finding the polygon index in set Q of blocking polys equal to the poly in which the principal node resides NOTE THIS IS PROBABLY NOT FUNCTIONAL RN 
                    poly_out = Q==k; %finding the polygon index in set K of blocking polys equal to the poly in which the terminal node resides
                    count = 1;
                    f = 1; % the number of blocking polygons after excluding the polygons in which the principal and terminal nodes reside
                    %blockersx = []; 
                    %blockersy = [];
                    polycost_Dijkstra = zeros(length(connection_blocked)-2, 1);
                    begindex_Dijkstra = zeros(length(connection_blocked)-2,1);
                    endex_Dijkstra = zeros(length(connection_blocked)-2,1);
                    blocking_poly_intersections = zeros(length(connection_blocked)-2,1);
                    while count <= length(connection_blocked)
                        if ismember(Q(connection_blocked(count)),Q(poly_in)) %if the node we're at is in the polygon
                            long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                            long2 = connections(i,2)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                            lat2 = connections(i,3)*pi/180;
                            diff_long = long2-long1; %finding the difference in longitude. 
                            diff_lat = lat2-lat1; %finding the difference in latitude.
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            if still_in_poly(connections(i,1),j) %if the node is still part of the original polygon, then include it
                                dist_total(i) = dist_total(i) + d*polycost(j); %adding this distance spent inside the home polygon to the total distance
                                weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(j);
                            else
                                dist_total(i) = dist_total(i) + d;
                                weighted_straight_line(i) = weighted_straight_line(i) + d;
                            end
                            
                            principal = [x_intersect(i,Q(connection_blocked(count)),1) y_intersect(i,Q(connection_blocked(count)),1)]; %store the outbound location as the boundary of the resident polygon
                            home_node = count; %store the index of the home node
                        elseif ismember(Q(connection_blocked(count)),Q(poly_out)) %if the node we're at is the terminal polygon, then we want to submit it as the boundary for Dijkstra
                            long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                            long2 = connections(i,5)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                            lat2 = connections(i,6)*pi/180;
                            diff_long = long2-long1; %finding the difference in longitude. 
                            diff_lat = lat2-lat1; %finding the difference in latitude.
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            if still_in_poly(connections(i,4),k) %if the terminal node is still in the polygon 
                                dist_total(i) = dist_total(i) + d*polycost(k); %adding this distance spent inside the home polygon to the total distance
                                weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(k); %add this distance spent inside the home polygon to the weighted straight line distance
                            else
                                dist_total(i) = dist_total(i) + d;
                                weighted_straight_line(i) = weighted_straight_line(i) + d;
                            end
                            terminal = [x_intersect(i,Q(connection_blocked(count)),1) y_intersect(i,Q(connection_blocked(count)),1)]; %store the final destination as the boundary of the resident node to polygon
                            away_node = count; %store the index of the away node
                        else
                            blockersx_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),1); nan]; %grab the second row
                            blockersy_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),2); nan]; %grab the third row
                            if f == 1 %if the first instance of discoverance of blocking polygon
                                blockersx = blockersx_append;
                                blockersy = blockersy_append;
                                begindex_Dijkstra(f) = 1;
                                endex_Dijkstra(f) = length(blockersx);
                            else %not the first instance
                                begindex_Dijkstra(f) = length(blockersx)+1;
                                blockersx = [blockersx; blockersx_append]; %piece together the vectors. 
                                blockersy = [blockersy; blockersy_append]; %piece together the vectors.
                                endex_Dijkstra(f) = length(blockersx);
                                %blockersx(f,:) = polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),2); %store for the shortest path algo the polygons that you'll be using, associated with the 
                                %blockersy(f,:) = polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),3);
                            end
                            polycost_Dijkstra(f) = polycost(Q(connection_blocked(count))); %store the associated polygon cost 
                            blocking_poly_intersections(f) = count; %store the count values so we can access easily the values of the polygon intersections. 
                            f = f+1;
                        end
                        count = count + 1;
                    end
                    count = 1; %counter for accessing blocking polygon intersections
                    f = 1; %counter for dealing with removing entries from Dijkstra since we don't want to tosic them completely because that skullfucks everything instead we use f as the queuer
                    while count <= length(blocking_poly_intersections) %now need to add up weighted straight line distances. How do I do this????
                        if length(blocking_poly_intersections) == 1
                            long1 = x_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            weighted_straight_line(i) = weighted_straight_line(i) + d; %add this straight line distance to the total distance
                            
                             
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1);
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2);
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1);
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2);
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            if lakejump %if we want to consider jumping this polygon
                                if d < resolution %if the distance is sufficiently small, we can ignore the weight of this in the straight weight calculation, and remove it from the shortest path blocking polygons.
                                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                                    blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %knock out all the points corresponding to this vector. 
                                    blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                                    polycost_Dijkstra(f) = [];
                                    begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1); % subtract off from the vector the length of the vector that was just removed. 
                                    endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                    begindex_Dijkstra(f) = []; %remove the begindices and endices because you don't need them anymore. 
                                    endex_Dijkstra(f) = [];
                                    f = f - 1; %because we straight up just removed a row we have to subtract one step so we don't skip over the next one. 
                                else
                                    weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add the weight of the blocking_poly_intersection index of the blocked_connections, and find that polygon in Q, then use polycost on that index
                                end
                            end
                            
                            %now add up the distance from the final
                            %blocking polygon to the away node. add 1 to
                            %count. 
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            long2 = x_intersect(i,Q(connection_blocked(away_node)),1)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            lat2 = y_intersect(i,Q(connection_blocked(away_node)),1)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            weighted_straight_line(i) = weighted_straight_line(i) + d;
                            
                            break
                        end
                        if count == 1 %we are between the boundary of the home polygon and the boundary of the first blocking polygon
                            long1 = x_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            weighted_straight_line(i) = weighted_straight_line(i) + d; %add this straight line distance to the total distance
                            
                            
                            %need to add up cost of first polygon as well
                            %as the tail after it. 
                            
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1);
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2);
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1);
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2);
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            if lakejump %if we want to consider jumping this polygon
                                if d < resolution %if the distance is sufficiently small, we can ignore the weight of this in the straight weight calculation, and remove it from the shortest path blocking polygons.
                                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                                    blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %knock out all the points corresponding to this vector. 
                                    blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                                    polycost_Dijkstra(f) = [];
                                    begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1); % subtract off from the vector the length of the vector that was just removed. 
                                    endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                    begindex_Dijkstra(f) = []; %remove the begindices and endices because you don't need them anymore. 
                                    endex_Dijkstra(f) = [];
                                    f = f - 1; %because we straight up just removed a row we have to subtract one step so we don't skip over the next one. 
                                else
                                    weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add the weight of the blocking_poly_intersection index of the blocked_connections, and find that polygon in Q, then use polycost on that index
                                end
                            end
                            
                            %now add the tail.... 
                                long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %end boundary of first blocking polygon
                                long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180; %beginning boundary of next blocking polygon 
                                lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                                lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180;
                                diff_long = long2 - long1;
                                diff_lat = lat2 - lat1;
                                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                                d = c * R; %step 3 of haversine function
                                weighted_straight_line(i) = weighted_straight_line(i) + d;
                            
                                
                        elseif count == length(blocking_poly_intersections) %we are between the boundary of the final blocking polygon and the boundary of the away polygon
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            long2 = x_intersect(i,Q(connection_blocked(away_node)),1)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            lat2 = y_intersect(i,Q(connection_blocked(away_node)),1)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            weighted_straight_line(i) = weighted_straight_line(i) + d;
                            
                            %need to add the cost of this final blocking
                            %polygon, if we don't lake jump it
                            
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1);
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2);
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1);
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2);
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            if lakejump %if we want to consider jumping this polygon
                                if d < resolution %if the distance is sufficiently small, we can ignore the weight of this in the straight weight calculation, and remove it from the shortest path blocking polygons.
                                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                                    blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %knock out all the points corresponding to this vector. 
                                    blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                                    polycost_Dijkstra(f) = [];
                                    begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1); % subtract off from the vector the length of the vector that was just removed. 
                                    endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                    begindex_Dijkstra(f) = []; %remove the begindices and endices because you don't need them anymore. 
                                    endex_Dijkstra(f) = [];
                                    f = f - 1; %because we straight up just removed a row we have to subtract one step so we don't skip over the next one. 
                                else
                                    weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add the weight of the blocking_poly_intersection index of the blocked_connections, and find that polygon in Q, then use polycost on that index
                                end
                            end
                        else %we are in betweeners. part 1 here is to add the distances between blocking polygons, part 2 is to add the distance inside the polygons. Hopefully polygons being within each other is not a problem because they will still determine internal distance because it's distances between each polygon. 
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %end boundary of first blocking polygon
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180; %beginning boundary of next blocking polygon 
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            weighted_straight_line(i) = weighted_straight_line(i) + d;
                            
                            
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180; %here is the distance inside the blocking polygon. Has an associated polycost
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %other boundary of the inside of the blocking polygon
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            if lakejump %if we want to consider jumping this polygon
                                if d < resolution %if the distance is sufficiently small, we can ignore the weight of this in the straight weight calculation, and remove it from the shortest path blocking polygons.
                                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                                    blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %knock out all the points corresponding to this vector. 
                                    blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                                    polycost_Dijkstra(f) = [];
                                    begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1); % subtract off from the vector the length of the vector that was just removed. 
                                    endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                    begindex_Dijkstra(f) = []; %remove the begindices and endices because you don't need them anymore. 
                                    endex_Dijkstra(f) = [];
                                    f = f - 1; %because we straight up just removed a row we have to subtract one step so we don't skip over the next one. 
                                else
                                    weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add the weight of the blocking_poly_intersection index of the blocked_connections, and find that polygon in Q, then use polycost on that index
                                end
                            end
                        end
                        count = count + 1;
                        f = f+1;
                    end
                    if shortestpath_bool %if the user wants to undergo a shortest path calculation
                        if dist_straight_line > shortestpath_lb && dist_straight_line < shortestpath_ub %if the straight line connection is within the bounds for analysis
                            if ~isempty(blockersx)
                                dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy,begindex_Dijkstra,endex_Dijkstra, polycost_Dijkstra,resolution,R);
                            else
                                dist_total(i) = weighted_straight_line(i);
                            end
                        else %otherwise the total distance is defined as the weighted straight line. 
                            dist_total(i) = weighted_straight_line(i);
                        end
                    else
                        dist_total(i) = weighted_straight_line(i); %if the user doesn't wish to do a shortest path algorithm, they can return simply the weighted straight line cost
                    end
                    %the nodes are in different polys but there is a blocking
                    %poly, then the distance is principal node ->
                    %boundary of poly1 -> boundary of poly2 -> other
                    %boundary of poly2 -> ... -> boundary of polyn ->
                    %terminal node. must decide whether it is cheaper to
                    %operate the line around the shortest boundary of the
                    %polygon or through the polygon based on polycost
                else 
                    %the nodes are in in different polys but there is no
                    %blocking poly,  then the distance is principal node ->
                    %boundary of poly1 -> boundary of poly2 -> terminal
                    %node SCENARIO 3
                    
                    [x_cross, y_cross] = polyxpoly(linspace(connections(i,2), connections(i,5),10),linspace(connections(i,3), connections(i,6),10),...
                        [polygons(begindex(j):endex(j),1); nan; polygons(begindex(k):endex(k),1)], [polygons(begindex(j):endex(j),2); nan; polygons(begindex(k):endex(k),2)]);
                    long1 = connections(i,2)*pi/180;
                    long2 = x_cross(1)*pi/180;
                    long3 = x_cross(2)*pi/180;
                    long4 = connections(i,5)*pi/180;
                    lat1 = connections(i,3)*pi/180;
                    lat2 = y_cross(1)*pi/180;
                    lat3 = y_cross(2)*pi/180;
                    lat4 = connections(i,6)*pi/180;
                    diff_long_1 = long1 - long2;
                    diff_long_2 = long2 - long3;
                    diff_long_3 = long3 - long4;
                    diff_lat_1 = lat1 - lat2;
                    diff_lat_2 = lat2 - lat3;
                    diff_lat_3 = lat3 - lat4;
                    a1 = sin(diff_lat_1/2)^2+cos(lat1)*cos(lat2)*sin(diff_long_1/2)^2; %step 1 of haversine function
                    c1 = 2*asin(min(1,sqrt(a1))); %step 2 of haversine function
                    d1 = c1 * R; %step 3 of haversine function
                    a2 = sin(diff_lat_2/2)^2+cos(lat2)*cos(lat3)*sin(diff_long_2/2)^2; %step 1 of haversine function
                    c2 = 2*asin(min(1,sqrt(a2))); %step 2 of haversine function
                    d2 = c2 * R; %step 3 of haversine function
                    a3 = sin(diff_lat_3/2)^2+cos(lat3)*cos(lat4)*sin(diff_long_3/2)^2; %step 1 of haversine function
                    c3 = 2*asin(min(1,sqrt(a3))); %step 2 of haversine function
                    d3 = c3 * R; %step 3 of haversine function
                    if ~still_in_poly(connections(i,1),j) %if we decided earlier that this node is no longer in this polygon we don't have to consider its cost in the final summation
                        if ~still_in_poly(connections(i,4),k) %if it's within the tolerance we can ignore the polycost
                            dist_total(i) = d1 + d2 + d3;
                            weighted_straight_line(i) = dist_total(i);
                        else %ignore d1 but not d3
                            dist_total(i) = d1 + d2 + d3*polycost(k);
                            weighted_straight_line(i) = dist_total(i);
                        end
                    else
                        if ~still_in_poly(connections(i,4),k) %ignore d3 but not d1
                            dist_total(i) = d1*polycost(j) + d2 + d3;
                            weighted_straight_line(i) = dist_total(i);
                        else
                            dist_total(i) = d1*polycost(j) + d2 + d3*polycost(k); %total distance is the sum of traveled time in each polygon and in free space with the associated penalty 
                            weighted_straight_line(i) = dist_total(i); %weighted_straight_line_distance is the same as total distance in this case. 
                        end
                    end
                    
                end
            elseif length(connection_blocked)>1
                %the terminal node is not in a poly, but there is a
                %blocking poly SCENARIO 4
                poly_in = Q==j; %finding the polygon index in set Q of blocking polys equal to the poly in which the principal node resides
                count = 1;
                f = 1;
                %blockersx = []; %creating the mini vector that is a subset containing the blocking polygons nicely lined up in a column vector
                %blockersy = [];
                polycost_Dijkstra = zeros(length(connection_blocked)-1, 1);
                blocking_poly_intersections = zeros(length(connection_blocked)-1,1);
                begindex_Dijkstra = polycost_Dijkstra;
                endex_Dijkstra = begindex_Dijkstra;
                
                while count <= length(connection_blocked) %for loop to sum together all of the polygon distances for this connection
                        if  ismember(Q(connection_blocked(count)),Q(poly_in)) %the principal node's polygon index in Q
                            long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                            long2 = connections(i,2)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                            lat2 = connections(i,3)*pi/180;
                            diff_long = long2-long1; %finding the difference in longitude. 
                            diff_lat = lat2-lat1; %finding the difference in latitude.
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            if ~still_in_poly(connections(i,1),j) %if it's not still in the polygon, we can ignore the weight
                                dist_total(i) = dist_total(i) + d;
                            else
                                dist_total(i) = dist_total(i) + d*polycost(j); %adding this distance spent inside the home polygon to the total distance
                                weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(j); %adding the weighted distance spent inside the home polygon to the total straight line weighted distance. 
                            end
                            principal = [x_intersect(i,Q(connection_blocked(count)),1) y_intersect(i,Q(connection_blocked(count)),1)]; %store the outbound location as the boundary of the resident polygon
                            home_node = count; %store the index of the home node
                        else % store the data for the rest of the polygons
                            blockersx_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),1); nan];
                            blockersy_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),2); nan];
                            if f == 1 %first blocking polygon instance
                                blockersx = blockersx_append;
                                blockersy = blockersy_append;
                                begindex_Dijkstra(f) = 1;
                                endex_Dijkstra(f) = length(blockersx);
                            else 
                                begindex_Dijkstra(f) = length(blockersx)+1;
                                blockersx = [blockersx; blockersx_append];
                                blockersy = [blockersy; blockersy_append];
                                endex_Dijkstra(f) = length(blockersx);
                            end
                            %blockersx(f,:) = polygonx(Q(connection_blocked(count)),:); %store for the shortest path algo the polygons that you'll be using, associated with the 
                            %blockersy(f,:) = polygony(Q(connection_blocked(count)),:);
                            polycost_Dijkstra(f) = polycost(Q(connection_blocked(count)));
                            blocking_poly_intersections(f) = count;
                            f = f+1;
                        end
                        count = count + 1;
                end
                terminal = [connections(i,5) connections(i,6)];
                count = 1;
                f = 1;
                while count <= length(blocking_poly_intersections)
                    if length(blocking_poly_intersections) == 1
                        long1 = x_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                        long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                        lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        weighted_straight_line(i) = weighted_straight_line(i) + d; %add this straight line distance to the total distance
                        
                        %add up distance inside polygon.........
                        long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        if lakejump %if we want to consider jumping this polygon
                            if d < resolution %if the distance is sufficiently small, we can ignore the weight of this in the straight weight calculation, and remove it from the shortest path blocking polygons.
                                weighted_straight_line(i) = weighted_straight_line(i) + d;
                                blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %knock out all the points corresponding to this vector. 
                                blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                                polycost_Dijkstra(f) = [];
                                begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1); % subtract off from the vector the length of the vector that was just removed. 
                                endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                begindex_Dijkstra(f) = []; %remove the begindices and endices because you don't need them anymore. 
                                endex_Dijkstra(f) = [];
                                f = f - 1; %because we straight up just removed a row we have to subtract one step so we don't skip over the next one. 
                            else
                                weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add the weight of the blocking_poly_intersection index of the blocked_connections, and find that polygon in Q, then use polycost on that index
                            end
                        end
                        
                        
                        %finally add the distance from the boundary of the
                        %terminal polygon to the terminal node
                         long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        long2 = connections(i,5)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        lat2 = connections(i,6)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        weighted_straight_line(i) = weighted_straight_line(i) + d; %distance between boundary of final bp and the terminal node. 
                        break
                    end
                        
                        
                        
                        
                        
                    if count == 1 %if we are at the first polygon then we can find the distance between the principal and first blocking polygon
                        long1 = x_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                        long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(home_node)),1)*pi/180;
                        lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        weighted_straight_line(i) = weighted_straight_line(i) + d; %add this straight line distance to the total distance
                        
                        %still need to calculate the distance inside the
                        %principal blocking polygon because it's the first
                        %iteration of the blocking polys 
                        long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        if lakejump %if we want to consider jumping this polygon
                            if d < resolution %if the distance is sufficiently small, we can ignore the weight of this in the straight weight calculation, and remove it from the shortest path blocking polygons.
                                weighted_straight_line(i) = weighted_straight_line(i) + d;
                                blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %knock out all the points corresponding to this vector. 
                                blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                                polycost_Dijkstra(f) = [];
                                begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1); % subtract off from the vector the length of the vector that was just removed. 
                                endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                begindex_Dijkstra(f) = []; %remove the begindices and endices because you don't need them anymore. 
                                endex_Dijkstra(f) = [];
                                f = f - 1; %because we straight up just removed a row we have to subtract one step so we don't skip over the next one. 
                            else
                                weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add the weight of the blocking_poly_intersection index of the blocked_connections, and find that polygon in Q, then use polycost on that index
                            end
                        end
                        
                        %now add the tail... 
                        if length(blocking_poly_intersections) > 1
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %end boundary of first blocking polygon
                            long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180; %beginning boundary of next blocking polygon 
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            weighted_straight_line(i) = weighted_straight_line(i) + d;
                        else
                            long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            long2 = connections(i,5)*pi/180;
                            lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                            lat2 = connections(i,6)*pi/180;
                            diff_long = long2 - long1;
                            diff_lat = lat2 - lat1;
                            a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                            c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                            d = c * R; %step 3 of haversine function
                            weighted_straight_line(i) = weighted_straight_line(i) + d;
                        end
                    elseif count == length(blocking_poly_intersections) %it's the distance between the final blocking polygon and the terminal node
                        long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        long2 = connections(i,5)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        lat2 = connections(i,6)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        weighted_straight_line(i) = weighted_straight_line(i) + d; %distance between boundary of final bp and the terminal node. 
                        %now have to calculate the distance inside the
                        %final blocking polygon
                        long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180; %introduce boundaries as within the shaperinos. 
                        long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        if lakejump %if we want to ignore the cost of this blocking polygon because it's smaller than the resolution
                            if d < resolution %if the distance is less than the resolution then we want to ignore it, both for the straight line calculation as well as Dijkstra's algorithm
                                weighted_straight_line(i) = weighted_straight_line(i) + d; %include only the straight line distance across this polygon.
                                blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %blockersx = setdiff(blockersx,blockersx(count,:)); % remove the blocking polygon from the blockers in the mesh since it will effectively be skipped over.
                                blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %blockersy = setdiff(blockersy,blockersy(count,:));
                                polycost_Dijkstra(f) = [];
                                begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                                endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                begindex_Dijkstra(f) = [];
                                endex_Dijkstra(f) = [];
                                f = f -1;
                            else
                                weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add this straight line distance to the total distance
                            end
                        end
                    else % we are in betweeners. need to add up the values in between blocking polygons as well as the time spent inside blocking polygons plus their polycost.
                        long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %the end of the first polygon along the line
                        long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180; %the beginning of the next polygon
                        lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        weighted_straight_line(i) = weighted_straight_line(i) + d; %add this straight line distance to the total distance
                        
                        long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180; %introduce boundaries as within the shaperinos. 
                        long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                        lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                        diff_long = long2 - long1;
                        diff_lat = lat2 - lat1;
                        a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                        c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                        d = c * R; %step 3 of haversine function
                        if lakejump %if we want to ignore the cost of this blocking polygon because it's smaller than the resolution
                            if d < resolution %if the distance is less than the resolution then we want to ignore it, both for the straight line calculation as well as Dijkstra's algorithm
                                weighted_straight_line(i) = weighted_straight_line(i) + d; %include only the straight line distance across this polygon.
                                blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %blockersx = setdiff(blockersx,blockersx(count,:)); % remove the blocking polygon from the blockers in the mesh since it will effectively be skipped over.
                                blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %blockersy = setdiff(blockersy,blockersy(count,:));
                                polycost_Dijkstra(f) = [];
                                begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                                endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                                begindex_Dijkstra(f) = [];
                                endex_Dijkstra(f) = [];
                                f = f - 1;
                            else
                                weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count)))); %add this straight line distance to the total distance
                            end
                        end
                    end
                    count = count + 1;
                    f = f+1;
                end
                if shortestpath_bool %if user wishes shortest path calculation 
                    if dist_straight_line(i) < shortestpath_ub && dist_straight_line(i) > shortestpath_lb %if this connection is within the range of desired checks
                        if ~isempty(blockersx)
                            dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy,begindex_Dijkstra, endex_Dijkstra, polycost_Dijkstra,resolution,R);
                        else
                            dist_total(i) = weighted_straight_line(i);
                        end
                    else
                        dist_total(i) = weighted_straight_line(i);
                    end
                else
                    dist_total(i) = weighted_straight_line(i);
                end
            else
                %the terminal node is not in a poly and there is not a
                %blocking poly SCENARIO 5
                x_cross = x_intersect(P(connection_blocked),Q(connection_blocked),1);
                y_cross = y_intersect(P(connection_blocked),Q(connection_blocked),1);
                long1 = connections(i,2)*pi/180;
                long2 = x_cross*pi/180;
                long3 = connections(i,5)*pi/180;
                lat1 = connections(i,3)*pi/180;
                lat2 = y_cross*pi/180;
                lat3 = connections(i,6)*pi/180;
                diff_long_1 = long2 - long1;
                diff_long_2 = long3 - long2;
                diff_lat_1 = lat2 - lat1; 
                diff_lat_2 = lat3 - lat2;
                a1 = sin(diff_lat_1/2)^2+cos(lat1)*cos(lat2)*sin(diff_long_1/2)^2; %step 1 of haversine function
                c1 = 2*asin(min(1,sqrt(a1))); %step 2 of haversine function
                d1 = c1 * R; %step 3 of haversine function
                a2 = sin(diff_lat_2/2)^2+cos(lat1)*cos(lat2)*sin(diff_long_2/2)^2; %step 1 of haversine function
                c2 = 2*asin(min(1,sqrt(a2))); %step 2 of haversine function
                d2 = c2*R; %step 3 of haversine function
                if ~still_in_poly(connections(i,1),j) %if it's within the tolerance then we ignore the penalty
                    dist_total(i) = d1 + d2;
                else
                    dist_total(i) = d1*polycost(j) + d2; %total distance is sum of penalized distance from principal node to the boundary, then straight line distance from there
                end
                weighted_straight_line(i) = dist_total(i);
            end
        
       
           
                
    elseif ismember(connections(i,4),I) %if the terminal node is in a polygon when the principal node is not
        k = J(find(I==connections(i,4),1)); %find the corresp polygon for the terminal node
        
        if length(connection_blocked) >= 2 %principal node not in poly, terminal node is, is there a blocking poly
            %SCENARIO 6
            poly_out = Q==k; %finding the polygon index in set Q of blocking polys equal to the poly in which the terminal node resides I DON'T THINK THIS FINNA WORK
            count = 1;
            f = 1;
            %blockersx = zeros(length(connection_blocked)-1, 4); %at worst it's 4 intersections per polygon if it passes through a corner or some shit. 
            %blockersy = zeros(length(connection_blocked)-1, 4);
            polycost_Dijkstra = zeros(length(connection_blocked)-1, 1);
            begindex_Dijkstra = polycost_Dijkstra;
            endex_Dijkstra = begindex_Dijkstra;
            blocking_poly_intersections = zeros(size(polycost_Dijkstra));
            while count <= length(connection_blocked) %over all the blocking polygons
                if ismember(Q(connection_blocked(count)), Q(poly_out)) %if the blocking poly is the away poly. checks for correspondence between the index of blocked connections and the index of the resident polygon from P and Q
                    
                    long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                    long2 = connections(i,5)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                    lat2 = connections(i,6)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    if still_in_poly(connections(i,4),k) %if the polygon is determined to still be in the polygon after the tolerance check
                        weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(k);
                    else
                        weighted_straight_line(i) = weighted_straight_line(i) + d;
                    end
                    terminal = [x_intersect(i,Q(connection_blocked(count)),1) y_intersect(i,Q(connection_blocked(count)),1)];
                    away_node = count;
                else
                    blockersx_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),1); nan];
                    blockersy_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),2); nan];
                    if f == 1 % if the first instance of blockers 
                        blockersx = blockersx_append;
                        blockersy = blockersy_append;
                        begindex_Dijkstra(f) = 1;
                        endex_Dijkstra(f) = length(blockersx);
                    else % it's not the first instance and we want to append them
                        begindex_Dijkstra(f) = length(blockersx)+1;
                        blockersx = [blockersx; blockersx_append];
                        blockersy = [blockersy; blockersy_append];
                        endex_Dijkstra(f) = length(blockersx);
                    end
                    %blockersx(f,:) = polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),1);
                    %blockersy(f,:) = polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),2);
                    polycost_Dijkstra(f) = polycost(Q(connection_blocked(count)));
                    blocking_poly_intersections(f) = count;
                    f = f+1;
                end
                count = count + 1;
            end
            principal = [connections(i,2) connections(i,3)];
            count = 1;
            f = 1;
            while count <= length(blocking_poly_intersections) %NEED TO CHANGE HERE THE INDICES FOR X INTERSECTIONS. DO THE SAME FOR SCENARIOS 2 AND 4. IN SCENARIO 4, LOOK FOR THE CASE WHERE WE ARE BETWEEN TEH FINAL BLOCKING POLY AND THE TERMINAL NODE AND ADD THAT TO SCENARIO 8. 
                if length(blocking_poly_intersections) == 1
                    %distance from principal node to only blocking poly. 
                    long1 = connections(i,2)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    lat1 = connections(i,3)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    weighted_straight_line(i) = weighted_straight_line(i)+d;
                    
                    
                    %distance inside the blocking polygon
                     long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180; %end boundary of first blocking polygon
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %beginning boundary of next blocking polygon 
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    if lakejump %if we want to jump it bc we it's smaller than the resolution
                        if d < resolution
                            weighted_straight_line(i) = weighted_straight_line(i) + d; %don't penalize the crossing
                            blockersx(begindex(f):endex(f)) = []; %get rid of the polygon as a blocker
                            blockersy(begindex(f):endex(f)) = [];
                            polycost_Dijkstra(f) = [];
                            begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                            endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                            begindex_Dijkstra(f) = [];
                            endex_Dijkstra(f) = [];
                            f = f - 1;
                        else
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count))));
                        end
                    end
                    
                    %finally go from boundary of final blocking polygon to
                    %boundary of terminal polygon
                    long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))), 2)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(away_node)), 1)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(away_node)),1)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                    break
                end
                    
                    
                    
                    
                    
                    
                if count == length(blocking_poly_intersections) %if we are at the final one, then its the final poly to the boundary of the away polygon also need to sum up the distance inside the final polygon. 
                    long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))), 2)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(away_node)), 1)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(away_node)),1)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                    
                    
                    %add up final amount 
                    long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180; %end boundary of first blocking polygon
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %beginning boundary of next blocking polygon 
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    if lakejump %if we want to jump it bc we it's smaller than the resolution
                        if d < resolution
                            weighted_straight_line(i) = weighted_straight_line(i) + d; %don't penalize the crossing
                            blockersx(begindex(f):endex(f)) = []; %get rid of the polygon as a blocker
                            blockersy(begindex(f):endex(f)) = [];
                            polycost_Dijkstra(f) = [];
                            begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                            endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                            begindex_Dijkstra(f) = [];
                            endex_Dijkstra(f) = [];
                            f = f - 1;
                        else
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count))));
                        end
                    end
                    
                    
                elseif count == 1 %we are at the principal node and need to connect to the first polygon. 
                    long1 = connections(i,2)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    lat1 = connections(i,3)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                    
                    %also need to perform the cross polygon check for this
                    %iteration of blocking polygons. If it jumps the lake,
                    %we can ignore it. 
                    long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    if lakejump %if we want to jump it bc we it's smaller than the resolution
                        if d < resolution
                            weighted_straight_line(i) = weighted_straight_line(i) + d; %don't penalize the crossing
                            blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %get rid of the polygon as a blocker
                            blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                            polycost_Dijkstra(f) = [];
                            begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                            endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                            begindex_Dijkstra(f) = [];
                            endex_Dijkstra(f) = [];
                            f = f - 1;
                        else
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count))));
                        end
                    end
                    
                    
                    %now need to add the tail... 
                    
                    long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180; %end boundary of first blocking polygon
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180; %beginning boundary of next blocking polygon 
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                    
                    
                else %we aren't at the final one and we need to sum both the distance from bp to bp as well as distance traveled internally in bp.
                    long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count+1))),1)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    weighted_straight_line(i) = weighted_straight_line(i) + d; %this calculation serves as inter blocking poly distances. 
                    
                    long1 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),1)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(blocking_poly_intersections(count))),2)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    if lakejump %if we want to jump it bc we it's smaller than the resolution
                        if d < resolution
                            weighted_straight_line(i) = weighted_straight_line(i) + d; %don't penalize the crossing
                            blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %get rid of the polygon as a blocker
                            blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                            polycost_Dijkstra(f) = [];
                            begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                            endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                            begindex_Dijkstra(f) = [];
                            endex_Dijkstra(f) = [];
                            f = f - 1;
                            
                        else
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(blocking_poly_intersections(count))));
                        end
                    end
                end
                count = count + 1;
                f = f + 1;
            end
            if shortestpath_bool
                if dist_straight_line(i) < shortestpath_ub && dist_straight_line(i) > shortestpath_lb
                    if ~isempty(blockersx)
                        dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy, begindex_Dijkstra, endex_Dijkstra, polycost_Dijkstra,resolution,R); %total distance is distance spent inside terminal polygon + shortest path algo from princip node to boundary of terminal node
                    end
                else
                    dist_total(i) = weighted_straight_line(i); %if we choose to ignore the shortest path calculation, then the total distance is equivalent to the weighted straight line
                end
            else
                dist_total(i) = weighted_straight_line(i);
            end
        else %principal node not in poly, terminal node is, but no blocking poly. then blocked is only the index of the crossage
            %SCENARIO 7
            x_cross = x_intersect(P(connection_blocked), Q(connection_blocked), 1);
            y_cross = y_intersect(P(connection_blocked), Q(connection_blocked), 1);
            long1 = connections(i,2)*pi/180;
            long2 = x_cross*pi/180;
            long3 = connections(i,5)*pi/180;
            lat1 = connections(i,3)*pi/180;
            lat2 = y_cross*pi/180;
            lat3 = connections(i,6)*pi/180;
            diff_long_1 = long2 - long1;
            diff_lat_1 = lat2 - lat1;
            diff_long_2 = long3 - long2;
            diff_lat_2 = lat3 - lat2;
            a1 = sin(diff_lat_1/2)^2+cos(lat1)*cos(lat2)*sin(diff_long_1/2)^2; %step 1 of haversine function
            c1 = 2*asin(min(1,sqrt(a1))); %step 2 of haversine function
            d1 = c1 * R; %step 3 of haversine function
            a2 = sin(diff_lat_2/2)^2+cos(lat2)*cos(lat3)*sin(diff_long_2/2)^2; %step 1 of haversine function
            c2 = 2*asin(min(1,sqrt(a2))); %step 2 of haversine function
            d2 = c2 * R; %step 3 of haversine function
            if ~still_in_poly(connections(i,4),k) %if we no longer consider it in poly, then the distance is not penalized
                dist_total(i) = dist_total(i) + d1 + d2;
                weighted_straight_line(i) = dist_total(i);
            else %we must consider the penalty
                dist_total(i) = dist_total(i) + d1 + d2*polycost(k);
                weighted_straight_line(i) = dist_total(i);
            end 
        end
    elseif ~isempty(connection_blocked) %neither node is in a polygon, but  there are blocking polygons
        %SCENARIO 8 
        principal = [connections(i,2) connections(i,3)];
        terminal = [connections(i,5) connections(i,6)];
        count = 1;
        
        %blockersx = zeros(length(connection_blocked), 4);
        %blockersy = zeros(length(connection_blocked), length(polygony));
        polycost_Dijkstra = zeros(length(connection_blocked), 1);
        begindex_Dijkstra = polycost_Dijkstra;
        endex_Dijkstra = begindex_Dijkstra;
        while count <= length(connection_blocked)
            blockersx_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),1); nan];
            blockersy_append = [polygons(begindex(Q(connection_blocked(count))):endex(Q(connection_blocked(count))),2); nan];
            if count == 1
                blockersx = blockersx_append;
                blockersy = blockersy_append;
                begindex_Dijkstra(count) = 1;
                endex_Dijkstra(count) = length(blockersx);
            else
                begindex_Dijkstra(count) = length(blockersx)+1;
                blockersx = [blockersx; blockersx_append];
                blockersy = [blockersy; blockersy_append];
                endex_Dijkstra(count) = length(blockersx);
                %blockersx(count,:) = polygonx(Q(connection_blocked(count)),:); %find the data points for the corresponding blocking polygon
                %blockersy(count,:) = polygony(Q(connection_blocked(count)),:);
                polycost_Dijkstra(count) = polycost(Q(connection_blocked(count)));
            end
            count = count + 1; 
        end
        count = 1;
        f = 1;
        while count <= length(connection_blocked)
            if length(connection_blocked) == 1
                long1 = connections(i,2)*pi/180;
                long2 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                lat1 = connections(i,3)*pi/180;
                lat2 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * R; %step 3 of haversine function
                weighted_straight_line(i) = weighted_straight_line(i) + d;
                
                
                %distance inside polygon
                long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                long2 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                lat2 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * R; %step 3 of haversine function
                if lakejump %if we want to jump it bc we it's smaller than the resolution
                    if d < resolution
                        weighted_straight_line(i) = weighted_straight_line(i) + d; %don't penalize the crossing
                        blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %get rid of the polygon as a blocker
                        blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                        polycost_Dijkstra(f) = [];
                        begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                        endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                        begindex_Dijkstra(f) = [];
                        endex_Dijkstra(f) = [];
                    else
                        weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(count)));
                    end
                end
                
                
                long1 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                long2 = connections(i,5)*pi/180;
                lat1 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                lat2 = connections(i,6)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d2 = c * R; %step 3 of haversine function
                weighted_straight_line(i) = weighted_straight_line(i) + d2;
                break
            end
                
                
                
                
                
                
            if count == 1 %we are at the first node so the distance is from the first node to the first polygon, then the distance inside the first blocking polygon. if lakejump, ignore the cost. 
                long1 = connections(i,2)*pi/180;
                long2 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                lat1 = connections(i,3)*pi/180;
                lat2 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * R; %step 3 of haversine function
                weighted_straight_line(i) = weighted_straight_line(i) + d;
                %also need to perform the cross polygon check for this
                %iteration of blocking polygons. If it jumps the lake,
                %we can ignore it. 
                long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                long2 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                lat2 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * R; %step 3 of haversine function
                if lakejump %if we want to jump it bc we it's smaller than the resolution
                    if d < resolution
                        weighted_straight_line(i) = weighted_straight_line(i) + d; %don't penalize the crossing
                        blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %get rid of the polygon as a blocker
                        blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                        polycost_Dijkstra(f) = [];
                        begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                        endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                        begindex_Dijkstra(f) = [];
                        endex_Dijkstra(f) = [];
                        f = f-1;
                    else
                        weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(count)));
                    end
                end
                
                %and now the tail...
                
                
                long1 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180; %end boundary of first blocking polygon
                long2 = x_intersect(i,Q(connection_blocked(count+1)),1)*pi/180; %beginning boundary of next blocking polygon 
                lat1 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                lat2 = y_intersect(i,Q(connection_blocked(count+1)),1)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d2 = c * R; %step 3 of haversine function
                weighted_straight_line(i) = weighted_straight_line(i) + d2;
                    
            elseif count == length(connection_blocked)%we want the distance between the final blocking poly and the terminal node
                    long1 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                    long2 = connections(i,5)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                    lat2 = connections(i,6)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    weighted_straight_line(i) = weighted_straight_line(i) + d;
                    
                    %add up distance traveled inside the final blocking
                    %polygon. 
                    long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                    long2 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                    lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                    lat2 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                    diff_long = long2 - long1;
                    diff_lat = lat2 - lat1;
                    a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                    c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                    d = c * R; %step 3 of haversine function
                    if lakejump %if we want to jump it bc we it's smaller than the resolution
                        if d < resolution
                            weighted_straight_line(i) = weighted_straight_line(i) + d; %don't penalize the crossing
                            blockersx(begindex_Dijkstra(f):endex_Dijkstra(f)) = []; %get rid of the polygon as a blocker
                            blockersy(begindex_Dijkstra(f):endex_Dijkstra(f)) = [];
                            polycost_Dijkstra(f) = [];
                            begindex_Dijkstra(f:end) = begindex_Dijkstra(f:end) - (endex_Dijkstra(f)-begindex_Dijkstra(f)+1);
                            endex_Dijkstra(f:end) = endex_Dijkstra(f:end) - (endex_Dijkstra(f) - begindex_Dijkstra(f)+1);
                            begindex_Dijkstra(f) = [];
                            endex_Dijkstra(f) = [];
                            f = f  - 1;
                            
                        else
                            weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(count)));
                        end
                    end
                    
                
            else %we're at an intermediate blocking polygon point :)
                long1 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                long2 = x_intersect(i,Q(connection_blocked(count+1)),1)*pi/180;
                lat1 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                lat2 = y_intersect(i,Q(connection_blocked(count+1)),1)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * R; %step 3 of haversine function
                weighted_straight_line(i) = weighted_straight_line(i) + d;
                
                long1 = x_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                long2 = x_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                lat1 = y_intersect(i,Q(connection_blocked(count)),1)*pi/180;
                lat2 = y_intersect(i,Q(connection_blocked(count)),2)*pi/180;
                diff_long = long2 - long1;
                diff_lat = lat2 - lat1;
                a = sin(diff_lat/2)^2+cos(lat1)*cos(lat2)*sin(diff_long/2)^2; %step 1 of haversine function
                c = 2*asin(min(1,sqrt(a))); %step 2 of haversine function
                d = c * R; %step 3 of haversine function
                if lakejump %if we wish to jump the lake
                    if d < resolution %if the distance is small enough we can ignore this polycost
                        weighted_straight_line(i) = weighted_straight_line(i) + d;
                    else
                        weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(count)));
                    end
                else
                    weighted_straight_line(i) = weighted_straight_line(i) + d*polycost(Q(connection_blocked(count)));
                end
            end
            count = count+1;
            f = f + 1;
        end
        if shortestpath_bool
            if dist_straight_line(i) < shortestpath_ub && dist_straight_line(i) > shortestpath_lb
                if ~isempty(blockersx)
                    dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy, begindex_Dijkstra, endex_Dijkstra, polycost_Dijkstra,resolution,R);
                else
                    dist_total(i) = weighted_straight_line(i);
                end
            else
                dist_total(i) = weighted_straight_line(i);
            end
        else
            dist_total(i) = weighted_straight_line(i);
        end
        
    else 
        %neither node is in a polygon and there are no intersecting polys
        %SCENARIO 9
        
        dist_total(i) = dist_straight_line(i); %total distance of connection is the same as the straight line distance
        weighted_straight_line(i) = dist_straight_line(i);
        
    end
    i = i + 1; %indexing!
end
        
                
 

















%for i = 1:length(I) %checking the num of connections for each node i in its corresp polygon
    %connection_indices(i) = find(node_from(:,1)==I(i));
 %   connections_from_size(i) = length(find(node_from(:,1)==I(i))); 
    
%end
%dist_node_in_poly = zeros(length(I),max(connections_from_size)); %preallocate matrix size based on i nodes in corresp polys, and the one with the most connections
%for i =1:length(I)%knowing which nodes are already in a polygon, find the distance of the connection inside the polygon
 %   node_de = find(node_from(:,1)==I(i)); %find index of node inside corresp polygon and all its connections
  %  node_a = node_to(node_de); %will return node indices of node i´s connections
   % node_de_lugar = [node_from(node_de(1),3),node_from(node_de(1),2)]; %go to first instance of node_from and collect those coordinates
   % node_a_lugar = [node_to(node_de,3),node_to(node_de,2)]; %collect coordinates of all corresponding connections to node i
   % for j = 1:length(node_de) %for all j connections of node i which resides in its ith polygon, need to check the intersections
    %    if isempty(setdiff(node_a(j),index_bus_in_poly(:,J(i)))) %if there are no elements in common between this node and the nodes in polygon i.e. this node is not in the same polygon
     %       [x_intersect,y_intersect] = polyxpoly([node_de_lugar(1),node_a_lugar(j,1)], ... will check for intersection point of polygon and connection
      %          [node_de_lugar(2), node_a_lugar(j,2)], polygonx(J(i,:)), polygony(J(i,:))); %checks polygon corresponding to node i 
       %     dist_node_in_poly = deg2km(sqrt((node_de_lugar(1)-x_intersect)^2+(node_de_lugar(2)-y_intersect)^2)); %node i in its corresp poly has j connections, whose dist in poly is recorded here
        %    connection_index = node_de(j); %returns index of the connection specified from node i to node j. works because node_de contains the indices within node_from, which has the same length as the total number of connections. the jth entry is the jth connection of type 'begins with node i', so it returns appropriately
        %
         %   dist_in_poly(connection_index,2) = dist_node_in_poly; %store distance penalty for final output
          %  dist_in_poly(connection_index,1) = J(i); %store polygon id for final output
          %  dist_in_poly(connection_index,4) = 1; %information re: number of nodes in the polygon in question
        %else
            %the polygon contains both nodes in the connection e.g. they
            %are both located in the same mountain range 
        %    dist_nodes_in_poly = deg2km(sqrt((node_de_lugar(1)-node_a_lugar(j,1))^2 + (node_de_lugar(2)-node_a_lugar(j,2))^2));
        %    connection_index = node_de(j);
         %   dist_in_poly(connection_index,2) = dist_nodes_in_poly;
         %   dist_in_poly(connection_index,1) = J(i);
         %   dist_in_poly(connection_index,4) = 2; %two nodes contained inside the polygon in question
       % end 
   % end
%end





%K = setdiff(all_nodes, I); %those nodes that are not in any polygon


%for i = 1:length(K)
 %   connections_K = find(node_from(:,1)==K(i));
  %  for k = 1:length(connections_K)
   %     for j = 1:length(row_polygon)
    %    [x_intersect,y_intersect] = polyxpoly([node_from_location(connections_K(k),2),node_to_location(connections_K(k),2)], ...
     %       [node_from_location(connections_K(k),1), node_to_location(connections_K(k),1)], polygonx(j,:),polygony(j,:)); % checks for intersections 
      %      if ~isempty(x_intersect) % if a connection exists
       %         x_cross_first = x_intersect(1);
        %        y_cross_first = y_intersect(1);
         %       x_cross_last = x_intersect(2);
          %      y_cross_last = y_intersect(2);
           %     dist_in_poly(connections_K(k),j) = deg2km(sqrt((x_cross_first - x_cross_last)^2 + (y_cross_first - y_cross_last)^2));
           % end
        %end
    %end
%end
            



%then establish node connections
%x_cross_first = zeros(length(node_from),length(polygonx));
%y_cross_first = zeros(length(node_from),length(polygonx));
%x_cross_last = zeros(length(node_from),length(polygonx));
%y_cross_last = zeros(length(node_from),length(polygonx));
%dist_in_poly = zeros(length(node_from),length(polygonx));
%index_line_in_poly = zeros(length(node_from),2);
%index_bus_in_poly = zeros(length(node_from),2);
%for i = 1:row_location
 %   for j = 1:row_polygon
  %      [x_intersect,y_intersect] = polyxpoly([node_from_location(i,2),node_to_location(i,2)], ...
   %         [node_from_location(i,1), node_to_location(i,1)], polygonx(j,:),polygony(j,:)); % checks for intersections 
    %    if ~isempty(x_intersect)
     %       x_cross_first(i,j) = x_intersect(1);
      %      y_cross_first(i,j) = y_intersect(1);
       %     if ~inpolygon(node_from_location(i,2),node_from_location(i,1),polygonx(j,:),polygony(j,:)) && ... %will check if both nodes are outside polygon
        %            ~inpolygon(node_to_location(i,2),node_to_location(i,1),polygonx(j,:),polygony(j,:))
         %       x_cross_last(i,j) = x_intersect(2);
          %      y_cross_last(i,j) = y_intersect(2);
           %     dist_in_poly(i,j) = deg2km(sqrt((x_cross_first - x_cross_last)^2 + (y_cross_first - y_cross_last)^2)); % distance that connection i spends in polygon j
            %    index_line_in_poly(i,:) = [i j]; %returns index of the connection i in polygon j 
            %%elseif index_bus_in_poly(i,j) %inpolygon(node_from_location(i,2),node_from_location(i,1),polygonx(j,:),polygony(j,:)) % if the node from is inpolygon
             %   dist_in_poly(i,j) = deg2km(sqrt((node_from_location(i,2) - x_intersect(1))^2 ...
       %             + (node_from_location(i,1) - y_intersect(1))^2)); 
        %        index_bus_in_poly(i,:) = [node_from(1,i),j];
         %       
          %  elseif ind%inpolygon(node_to_location(i,2),node_to_location(i,1),polygonx(j,:),polygony(j,:)) %if the node to is inpolygon
           %     dist_in_poly(i,j) = deg2km(sqrt((node_to_location(i,2) - x_intersect(1))^2 ...
            %        + (node_to_location(i,1) - y_intersect(1))^2));
             %   index_bus_in_poly(i,:) = [node_to(1,i),j];

            %else %both nodes are in polygon
            %    dist_in_poly(i,j) = deg2km(sqrt((node_from_location(i,2) - node_to_location(i,2))^2 + ...
                %    (node_from_location(i,1) - node_to_location(i,1))^2));
 %               index_bus_in_poly(i,:) = [1e+6,j]; %using 1e+6 to signify that both nodes are in the polygon
  %              
   %         end
    %    end
    %end
%end
%dist_in_poly = sparse(dist_in_poly);
%index_bus_in_poly = sparse(index_bus_in_poly);
%index_line_in_poly = sparse(index_line_in_poly);
toc
end

