/* eslint-disable */

import { Box, Button, Typography } from "@mui/material"
import { useEffect, useRef, useState } from "react"
import { FaLessThanEqual } from "react-icons/fa"
import { useHistory } from "react-router-dom"
import { colors } from "shared/theme/colors"
import CustomButton from "components/CustomButton";

// Outlined Button specific to Models and Atlases
/**
 * OutlinedButton
 * Specific to Model and Atlas Cards
 * @param content text content to be displayed
 * @param link button href link  
 */
export const OutlinedButton = ({ content, link = null, onClick, bg, color, bgHover, colorHover }) => {
  return (
    <Button
      variant="outlined"
      disableRipple
      href={link ? `#${link}` : null}
      onClick={onClick}
      sx={{
        color: color ? color : colors.primary[100], borderWidth: "2px", borderColor: bg ? bg : colors.primary[100], borderRadius: "1.2rem", justifyContent: 'center', textAlign: 'center',
        ":hover": { borderColor: colorHover ? colorHover : '#01579B', borderWidth: "2px", backgroundColor: bgHover ? bgHover : colors.primary[100], color: colorHover ? colorHover : '#01579B', transition: '0.5s' }
      }}
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
export const ModelCard = ({ width = "100%", height = "100%", title, description, onClick, learnMoreLink, disabled = false, onSelect, selected=true }) => {

  const [hover, setHover] = useState(false)
  const ref = useRef()
  const [flexDir, setFlexDir] = useState("column")

  useEffect(() => {

    // checks if the parent element is wider than it is longer
    // if it is, converts the flex direction to row
    if (ref.current.clientWidth > ref.current.clientHeight) {
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
          cursor: disabled ? "pointer": "initial",
          display: "flex",
          flexDirection: "column",
          justifyContent: "center"
        }}
        onMouseEnter={() => setHover(true)}
        onMouseLeave={() => setHover(false)}
      >
        {
          !disabled && hover &&
          <Box
            style={{ background: "linear-gradient(#4F83CC, #01579B)" }}
            sx={{
              position: "absolute",
              width: "100%",
              height: "100%",
              borderRadius: "1.2rem",
              p: "1rem",
              opacity: 0.9,
              boxShadow: "0px 4px 6px 0px rgba(1, 87, 155, .20), 0px 0px 1px 0px rgba(1, 87, 155, .32)"
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
              <OutlinedButton content="Select" onClick={onSelect} />
              <OutlinedButton content="Learn More" link={learnMoreLink} onClick={(e) => e.stopPropagation()}/>
              {
                disabled &&
                <Typography sx={{ color: colors.primary[900], fontSize: "12px", textDecoration: 'underline', textAlign: "center" }}>
                  Model unavailable for mapping
                </Typography>
              }
            </Box>
          </Box>
        }
        {/* DISABLED OVERLAY BOX */}
        {
          disabled &&
          <Box
            style={{ background: "linear-gradient(#e7e7e7, #d0d0d0)" }}
            sx={{
              position: "absolute",
              width: "100%",
              height: "100%",
              borderRadius: "1.2rem",
              p: "1rem",
              opacity: 0.9,
            }}>
            <Typography sx={{ position: "absolute", fontSize: "12px", fontWeight: "bold", color: colors.neutral[900], textAlign: "center", left: "28%" }}>
              Not Compatible
            </Typography>
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

              { hover && <OutlinedButton content="Learn More" link={learnMoreLink} bg={colors.neutral[800]} color={colors.neutral[800]} onClick={(e) => e.stopPropagation()} bgHover={colors.neutral[100]} colorHover={colors.neutral[600]}/> }
              {/* <CustomButton type="outlined" href={learnMoreLink ? `#${learnMoreLink}` : null} onClick={(e) => e.stopPropagation()}>Learn More</CustomButton> */}
            </Box>
          </Box>
        }
        { 
          !disabled &&
          <Box sx={{
          width: '100%',
          height: "100%",
          display: "flex",
          flexDirection: "column",
          p: "1.2rem",
          boxShadow: hover ? "none" : "0px 0px 10px rgba(0, 0, 0, 0.15)",
          borderRadius: "1.2rem",
          borderStyle:"solid",
          borderColor: selected ? "#008BF5" : 'transparent',
          borderWidth:"4px",
        }}
        >
          <Typography sx={{ fontSize: "1.4rem", fontWeight: "bold" }}>{`Model ${title}`}</Typography>
          <Typography sx={{ fontSize: "1rem", color: colors.neutral[800] }}>{description}</Typography>
        </Box>}
        { 
          disabled &&
          <Box sx={{
            width: '100%',
            height: "100%",
            display: "flex",
            flexDirection: "column",
            p: "1.2rem",
            borderRadius: "1.2rem",
          }}
          >
            <Typography sx={{ fontSize: "1.4rem", fontWeight: "bold" }}>{title}</Typography>
            <Typography sx={{ fontSize: "1rem", color: colors.neutral[800] }}>{description}</Typography>
          </Box>
        }
      </Box>
    </Box>
  )
}
