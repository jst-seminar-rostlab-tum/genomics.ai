import { Box, Typography } from '@mui/material'
import { useEffect, useRef, useState } from 'react'
import { colors } from 'shared/theme/colors'

export default function DashboardHeader(props){

    const {title, fontSize, children, width="100%", height="100%"} = props

    const titleRef = useRef()
    const totalRef = useRef()
    const [titleWidth, setTitleWidth] = useState(0)

    useEffect(()=>{
        setTitleWidth(totalRef.current.clientWidth - titleRef.current.clientWidth)
    }, [])

    return (
        <Box ref={totalRef} sx={{width, height, position: "relative", display: "flex", flexDirection: "row", alignItems: "center"}}>

            <Box ref={titleRef}><Typography fontWeight="bold" fontSize={fontSize ? fontSize : "2.2em"}>{title}</Typography></Box>
            <Box sx={{width: titleWidth, height: "100%"}}>{children}</Box>

            <Box sx={{position: "absolute", height: "2px", width: "100%", top: "100%", left: "0%", border: `1px solid ${colors.neutral[500]}`}} />
        </Box>
    )
}