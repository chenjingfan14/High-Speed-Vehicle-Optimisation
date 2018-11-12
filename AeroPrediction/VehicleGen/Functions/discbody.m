function new = discbody(old)

new = old.Points;
new.name = "aftbody";

target = 1.1;

xPanels = round(new.x(end)/target);

diff = new.x(2)-new.x(1);
new.x = (new.x(1):diff/xPanels:new.x(2))';

new = xyztopoints(new);