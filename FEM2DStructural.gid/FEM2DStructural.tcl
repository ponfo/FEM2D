proc InitGIDProject { dir } {
    Splash $dir
}
proc Splash { dir } {
    global GIDDEFAULT
    set VersionNumber "FEM2D - v0.9 - 2019"
    if { [.gid.central.s disable windows] } { return }
    if { [ winfo exist .splash]} {
	destroy .splash
	update
    }
    toplevel .splash        
    set im [image create photo -file [ file join $dir FEM2D.jpg]]
    set x [expr [winfo screenwidth .splash]/2-[image width $im]/2]
    set y [expr [winfo screenheight .splash]/2-[image height $im]/2]

    wm geom .splash +$x+$y
    wm transient .splash .gid
    wm overrideredirect .splash 1
    pack [label .splash.l -image $im -relief ridge -bd 2]    
    label .splash.lv -text $VersionNumber -font "times 7"
    place .splash.lv -x 500 -y 253
    bind .splash <1> "destroy .splash"
    bind .splash <KeyPress> "destroy .splash"
    raise .splash .gid
    grab .splash
    focus .splash
    update
    after 7000 "if { [ winfo exist .splash] } { 
	destroy .splash
    }"
}

