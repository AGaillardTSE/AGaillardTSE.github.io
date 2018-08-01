---
layout: indexing
---

# This is a test for index

<script>

var pullfiles=function(){ 
    // love the query selector
    var fileInput = document.querySelector("#myfiles");
    var files = fileInput.files;
    // cache files.length 
    var fl = files.length;
    var i = 0;

    while ( i < fl) {
        // localize file var in the loop
        var file = files[i];
        alert(file.name);
        i++;
    }    
}

// set the input element onchange to call pullfiles
document.querySelector("#myfiles").onchange=pullfiles;

//a.t
</script>

