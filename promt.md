okay, it runs. 

the problem is that the bridge creation and also the hole creation in the slabs is not correct. 

i now try to explain you in great detail what i want ans why i want it. listen and think strong. It is complicated.


goal: creating a mesh, that when printed in lw pla, a lightweight internal structure is printed without generating the internal structure in the slicer. The reason why we do that is that we can control the the internal structure such that the part is printed as close to vase mode as possible. we want to do this, because lw pla expands and cannot be retracted without stringing. 

How we can create the internal structure: immagine you have a really simple mesh. now you cut out a a slice from the border of the solid to approximatly the middle of the solid all the way from -z to +z of the solid. what happens now when we print the solid in the z direction, is that the outer wall of the solid is printed, including the thin slice we substracted from the solid. And like that we have created an internal structure that can be printed in lw pla. now lets go a step further. We do not want internal structures that go just to the center of the solid, because they do not really improve strength. instead wemake the structure throughout the solid: we substarct from the same solid a slice from -z to z but this this it goad through the whole soid. we basically have two parts, that cannot be printed without retraction. this is where the bridges come into play. the small bridges prevent the solid to be cut into multiple parts from the slab and provide a path for the printer to go. 

then we also make holes in the support structure such that the internal structure is lighter. 

practically we create a slab mesh, a bridgee mesh and hole mesh. 

from the slab mesh we can substract the hole mesh and the bridge meash. this resulting mesh is then substracted from the solid. and now we have waht we want. 


I will now explain how geometrically we want to create the bridge mesh:
Immagine that each slab needs to be cut into two along the z direction. we slice the solid and the slabs in the z-direction in x steps with dz. this gives us a whire for each slab and the solid. now we itterate over each slabwhire (lets represent the slab as a line). the slab whire has a certain length that is defined by the intersection eith the solid. we calculate the middle of this linesegment and this is where the bridge must go throug. we collect all those points for each slab in the z direction. so for each slab we have then a set of points, which can be used to define a bridge. the points of one bridge for one slab can be sonnected to a whire. the whire can be extruded to a surface the perpendicularly intersects the two surfaces of the slab. then we can extrude it into the other direction to give it the width. 

I will now explain how geometrically we want to create the hole mesh:
for this we need slab sections. the slab sections we can get by intersection cutting the slabs with each other. each slab section gets a hole. heach hole is defined by four points.
first we try to find the lowest and the highest  point(s) in the z direction of each slab section. if there is only one point at the highes /lowest z. those points are used for the hole definition. if there are multiple points (they are on a straigth line segment), the middle of the line segment are the points that define the hole. now we have two points of the four. the other thwo points are gotten by connecting the thwo points we already have and finding the modpoint. this midpoint can be extruded perpendigular to the line section within the slab section. the intersections of this line with the slabsection gives us the other two point. those four points define a surface. that can be ectruded to the slabwidth. and this is the hole. 


lets try to go stpe by step. i suggest we firrst create the slabs, the bridges and holes, only as surfaces. we can then gut the slab surfaces with the hole surface snd the bridge surfaces. this gives us the internal structure which we can extrude, to give us a thickness. 

we now simplyfy first the code such that it only creates the slabs as surfaces, that are already cut to the wing. 



what this gode is supposed to do:
create an internal structure (ribs) of the solid as a mesh. 
i have two sets of ribs that intersect each other. 
I need ne cut the ribs with each other and get the parts of the ribs. 



def create_bridges_analytical()
    bridges = []
    for rib in ribs:
        points_for_bridge = []
        rib_surface = build_analytical_ribsurface_equation
        z_min, z_max = get_bound(rib) # the distance only in the z direction
        for z in range(z_min, z_max, dz):
            slice = slice_slice_wing_at(z)
            rib_intersection = get_rib_intersection_in_z(z, ribsurface, slice)
            start_point, end_point = get_start_and_endpoint_of_line_segemnt(rib_intersection) # here you did not understand
            center_point = get_center_from_two_points(strt_point, end_point)
            point1 = offset_point_on_line(centerpoint, rib_intersection, 0.5*rib_height)
            point2 = offset_point_on_line(centerpoint, rib_intersection, -0.5*rib_height)
            points_for_bridge.append([p1, p2])
        bridge = make bridge_mesh_from_points(points_for_bridge)
        bridges.append(bridge)


def get_start_and_endpoint_of_line_segemnt(rib_intersection)
    # rib intersection is a line segment representing the rib intersection with the wing at the z slice.
    # the start and endpoint is the sart and endpoint of the ribintersection with the wing and not the start and end of the slab itself.
    # so the start and endpoint lie on the wingsurface.
    return start_point, end_point





def create_holes_analytical(point_condition:type funciton)
    bridges = []
    for rib_part in rib_parts:
        points_for_bridge = []
        rib_surface = build_analytical_ribsurface_equation()
        z_min, z_max = get_bound(rib_part)
        rib_height_in_z = z_max - z_min
        hole_height_in_z = z_max - z_min - 2 * hole_marging
        for z in range(z_min, z_max, dz):
            if z > hole_marging and z < (rib_height_in_z-hole_margin):
                slice = slice_slice_wing_at(z)
                rib_intersection = get_rib_intersection_in_z(z, ribsurface, slice)
                start_point, end_point = get_start_and_endpoint_of_line_segemnt(rib_intersection)
                center_point = get_center_from_two_points(strt_point, end_point, end_point, hole_marging )
                point1, point2 = make_points_from_condition(point_condition, rib_intersection, start_point, center_point, hole_height_in_z)
                points_for_bridge.append([p1, p2])
            
        bridge = make bridge_mesh_from_points(points_for_bridge)
        bridges.append(bridge)


def make_points_from_condition(point_condition, rib_intersection, start_point, center_point, rib_height_in_z, hole_marging, hole_height_in_z):
    distance_between_start_end = get_dist(start_point, center_point) - 2 * hole_marging
    if distance_between_start_end < 0:
        distance_between_start_end = 0
    x = (rib_height_in_z-hole_marging) / hole_height_in_z # this is a relative measure between 0 and 1 on how far we are in the rib slicing
    dist_relative = point_condition(x=x)
    point1 = center_point + (distance_between_start_end/2)*dist_relative
    point2 = center_point - (distance_between_start_end/2)*dist_relative
    return point1, point2

def point_condition(x):
    """x and y must be numbers between 0 and 1"""
    y = 0.5
    return y