import { useLocation } from 'react-router-dom';

import { Box, Typography } from '@mui/material'
import ArrowForwardIosIcon from '@mui/icons-material/ArrowForwardIos';
import { colors } from 'shared/theme/colors'
import { Link } from 'react-router-dom'

/**
 * Breadcrumb needs two parameters:
 * 1. fontsize
 * 2. path
 *     Format of Path should be like:
 *     - /explore/atlases
 */
export default function Breadcrumb({ fontSize, actions }) {

    let position = useLocation()

    let elements = position.pathname.split('/')

    const iconSize = fontSize * 0.7

    const rebuildLink = (index) => {
        let accmulator = ""
        for(let i = 1; i<=index; i++) accmulator = `${accmulator}/${elements[i]}`
        return accmulator
    }

    const handleOnClick = (element) => {
        if(!actions) return;
        const action=actions[element]
        if(action && typeof(action)==="function") action()
    }
    
    const generate = (element, index) => index === 0 ? <div key={index} /> : (
        <Box sx={{display: "flex", flexDirection: "row", alignItems: "center"}} key={index}>
            { index != elements.length - 1 
                ? <Box component={Link} to={rebuildLink(index)} sx={{cursor: "pointer", textDecoration: "none"}} onClick={()=>handleOnClick(element)}><Typography fontWeight="bold" fontSize={`${fontSize}em`} sx={{color: colors.neutral[500], textTransform: "capitalize"}}>{element}</Typography></Box> 
                : <Box><Typography fontWeight="bold" fontSize={`${fontSize}em`} sx={{color: colors.primary[800], textTransform: "capitalize"}}>{element}</Typography></Box>}
            { index === elements.length - 1 ? <></> : <ArrowForwardIosIcon sx={{width: `${fontSize*iconSize}em`, height: `${fontSize*iconSize}em`, margin: `0em ${fontSize*0.2}em`}}/> }
        </Box>
    )

    return (
        <Box sx={{display: "flex", flexDirection: "row", alignItems: "center"}}>
            {
                elements.map((element, index) => generate(element, index))
            }
        </Box>
    )
}