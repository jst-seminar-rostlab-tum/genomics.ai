/* eslint-disable */

import { Box, Button, Typography } from "@mui/material"
import { useEffect, useRef, useState } from "react"
import { colors } from "shared/theme/colors"

// Outlined Button specific to Models and Atlases
/**
 * OutlinedButton
 * Specific to Model and Atlas Cards
 * @param content text content to be displayed
 * @param link button href link  
 */
export const OutlinedButton = ({ content, link }) => {
 return (
   <Button
    variant="outlined"
    disableRipple
    href={`/${link}`}
    sx={{ color: "white", borderWidth: "2px", borderColor: "white", borderRadius: "1.2rem", ":hover": { borderColor: colors.secondary1[500], borderWidth: "2px" } }}
   >
    {content}
   </Button>
 )
}

/**
 * ModelCard
 * @param width default value 100% of parent
 * @param height default value 100% of parent
 * @param title
 * @param description
 * @param mapLink onHover button Map url
 * @param learnMoreLink onHover button Learn More url
 */
export const ModelCard = ({ width = "100%", height = "100%", title, description, mapLink, learnMoreLink}) => {

  const [hover, setHover] = useState(false)
  const ref = useRef()
  const [flexDir, setFlexDir] = useState("column")

  useEffect(() => {

    // checks if the parent element is wider than it is longer
    // if it is, converts the flex direction to row
    if(ref.current.clientWidth > ref.current.clientHeight) {
      setFlexDir("row")
    }
  }, [])

  return (
    <Box
      sx={{
        width: width,
        height: height
      }}
    >
      <Box
        ref={ref}
        sx={{ 
          width: "100%", 
          height: "100%", 
          position: "relative", 
          cursor: "pointer", 
          display: "flex", 
          flexDirection: "column", 
          justifyContent: "center"
        }}
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
        >  
        {
          hover && 
          <Box 
          style={{ background: "linear-gradient(#4F83CC, #01579B)" }} 
          sx={{ 
            position: "absolute", 
            width: "100%", 
            height: "100%", 
            borderRadius: "1.2rem",
            p:"1rem",
            opacity: 0.9
          }}>
            <Box
              sx={{
                display: "flex",
                flexDirection: flexDir,
                gap: "1rem",
                // styles below center the buttons to the very center of the component
                position: "absolute",
                top: "50%",
                left: "50%",
                transform: "translate(-50%, -50%)"
              }}  
              >
              <OutlinedButton content="Map" link={mapLink}/>
              <OutlinedButton content="Learn More" link={learnMoreLink}/>
            </Box>
          </Box>
        }
      <Box sx={{ 
        width: '100%',
        height: "100%",
        display: "flex", 
        flexDirection: "column", 
        p: "1.2rem", 
        boxShadow: "0px 0px 10px rgba(0, 0, 0, 0.15)", 
        borderRadius: "1.2rem",
      }}
      >
        <Typography sx={{ fontSize: "1.4rem", fontWeight: "bold" }}>{title}</Typography>
        <Typography sx={{ fontSize: "1rem", color: colors.neutral[800] }}>{description}</Typography>
      </Box>
      </Box>
    </Box>
  )
}