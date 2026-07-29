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