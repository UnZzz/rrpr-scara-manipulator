@tool
extends Node3D

@onready
var x_label : Label3D = $Node3D/Label3D

@onready
var y_label : Label3D = $Node3D/Label3D2

@onready
var z_label : Label3D = $Node3D/Label3D3

@export
var i : String


func _process(delta: float) -> void:
	x_label.text = "X" + i
	y_label.text = "Y" + i
	z_label.text = "Z" + i
