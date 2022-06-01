import { useLocation, Link } from 'react-router-dom';

import { Box, Typography } from '@mui/material';
import ArrowForwardIosIcon from '@mui/icons-material/ArrowForwardIos';
import { colors } from 'shared/theme/colors';

/**
 * Breadcrumb needs the parameter fontsize
 * 
 * It may also accept another parameter which is actions
 * 
 * the actions is the object, which summarize the action when clicking the corresponding element
 * 
 * for example:
 * 
 * the path is /explore/atlases
 * 
 * if you want to do something additionally than just redirect
 * 
 * you should give an object like this:
 * 
 * {
 *      explore: ()=>{something you want to do}
 * }
 * 
 */
export default function Breadcrumb({ fontSize, actions, path, elems, p = 0 }) {

  // let position = useLocation()

  let elements = elems ? elems : path.split('/')

  const iconSize = fontSize * 0.7

  const rebuildLink = (index) => {
    let accmulator = ""
    for (let i = 1; i <= index; i++) accmulator = `${accmulator}/${elements[i]}`
    return accmulator
  }

  const handleOnClick = (element) => {
    if (!actions) return;
    const action = actions[element]
    if (action && typeof (action) === "function") action()
  }

  const generate = (element, index) => index === 0 ? <div key={index} /> : (
    <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center" }} key={index}>
      {index != elements.length - 1
        ? <Box component={Link} to={rebuildLink(index)} sx={{ cursor: "pointer", textDecoration: "none" }} onClick={() => handleOnClick(element)}><Typography fontWeight="bold" fontSize={`${fontSize}em`} sx={{ color: colors.neutral[500], textTransform: "capitalize" }}>{element}</Typography></Box>
        : <Box><Typography fontWeight="bold" fontSize={`${fontSize}em`} sx={{ color: colors.primary[800], textTransform: "capitalize" }}>{element}</Typography></Box>}
      {index === elements.length - 1 ? <></> : <ArrowForwardIosIcon sx={{ width: `${fontSize * iconSize}em`, height: `${fontSize * iconSize}em`, margin: `0em ${fontSize * 0.2}em` }} />}
    </Box>
  )

  return (
    <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", paddingLeft: p, width: "100%" }}>
      {
        elements.map((element, index) => generate(element, index))
      }
    </Box>
  )
}
