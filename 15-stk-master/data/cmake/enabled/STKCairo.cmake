if(${STK_OPT_CAIRO_DISABLE})
	set(STK_CAIRO_ENABLED 0)
else(${STK_OPT_CAIRO_DISABLE})
	add_definitions(-DCAIRO_ENABLED)
	set(STK_CAIRO_ENABLED 1)
endif(${STK_OPT_CAIRO_DISABLE})