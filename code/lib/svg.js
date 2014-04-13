function generate_xy(num, maxx, maxy) {
    var xys = [];
    for(var i=0; i<num; i++) {
        xys.push({x: Math.random()*(maxx-20)+10, y: Math.random()*(maxy-20)+10});
    }
    return xys;
}

function connect_circles(circles) {
    var lines = [];
    circles.forEach( function(a) {
        circles.forEach( function(b) {
            if(a.x == b.x && a.y == b.y) return;
            lines.push({x1: a.x, y1: a.y, x2: b.x, y2: b.y});
        })
    })
    return lines;
}
