## Extract the row numbers for doubly imaged sources
function abs(x)
{ 
    return x<0 ? -x:x
}
{
    if($NF~/images/ && $1==2) 
    {
	bool=1;
	next;
    }; 
    
    if(bool==1 && $NF~/parity/ && abs($3)>1.0) 
    {
	print NR-1; 
	bool=0;
    } 
    else
    {
	bool=0;
    }
}
