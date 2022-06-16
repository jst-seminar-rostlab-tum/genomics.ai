import { Box, Typography } from '@mui/material'
import { useEffect, useRef, useState } from 'react'
import { colors } from 'shared/theme/colors'

/**
 * Dashboard Header
 * 
 * the size and font size and the title can be customized
 * 
 * I leave the area right next to the title free and If you want to add something there
 * you can use children, the idea of this is that the additional elements like add button 
 * and the textfield may have lots of interaction with parent component and it's better to let user create it
 * 
 * example:
 * <DashboardHeader title="Title" fontSize="3.2em" width="800px" height="150px">
 * 
 * <Box sx={{position: "relative", top: "13px", left: "20px", display: "flex", flexDirection: "row", alignItems: "center"}}>
 *     <Typography fontSize="1em">Technische Universität München</Typography>
 * 
 *     <Box sx={{position: "relative", width: "100px", height: "30px", marginLeft: "10px"}}>
 *         <Tag content="Public" variant="primary-default" />
 *     </Box>
 * </Box>
 * 
 * <Button sx={{marginLeft: "220px"}}>Add</Button>
 * </DashboardHeader>
 * 
 */
export default function DashboardHeader(props){

    const {title, fontSize, children, width="100%", height="100%"} = props

    const titleRef = useRef()
    const [titleWidth, setTitleWidth] = useState(0)

    useEffect(()=>{
        setTitleWidth(titleRef.current.clientWidth)
    })

    return (
        <Box sx={{width, height, position: "relative", display: "flex", flexDirection: "row", alignItems: "center"}}>

            <Box ref={titleRef}><Typography fontWeight="bold" fontSize={fontSize ? fontSize : "2.2em"}>{title}</Typography></Box>
            <Box sx={{position: "relative", width: `calc(100% - ${titleWidth}px)`, height: "100%", display: "flex", flexDirection: "row", alignItems: "center"}}>{children}</Box>

            <Box sx={{position: "absolute", height: "2px", width: "100%", top: "100%", left: "0%", border: `1px solid ${colors.neutral[500]}`}} />
        </Box>
    )
}