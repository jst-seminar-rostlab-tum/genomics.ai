/* eslint-disable */

import { Box, IconButton, Typography } from "@mui/material"
import CopyrightIcon from '@mui/icons-material/Copyright';
import { Link } from "react-router-dom"

import { colors } from "shared/theme/colors";
import logo from 'assets/logo.svg';

function Footbar(props) {
  return (
    <Box {...props} sx={{
      margin: "auto",
      display: "flex",
      flexDirection: "row",
      flexWrap: "warp",
      justifyContent: "space-between",
      alignItems: "center"
    }}
    ></Box>
  )
}

function Navlink(props) {
  return (
    <Box {...props}>
      <Typography {...props} sx={{
        fontWeight: "bold",
        fontSize: { xs: "0.6em", sm: "0.8em", md: "0.8em", lg: "1.2em", xl: "1.2em" },
        color: "black",
        ":hover": {
          color: colors.primary[400]
        }
      }}
      ></Typography>
    </Box>
  )
}

function LinkBox(props) {
  return (
    <Box {...props} component={Link} style={{ textDecoration: "none" }}></Box>
  )
}

export default function Footer({ sx }) {
  return (
    <Box
      sx={{
        position: "relative",
        bgcolor: colors.neutral[200],
        width: "100%",
        height: "8vh",
        display: "flex",
        flexDirection: "row",
        alignItems: "center",
        justifyContent: "center",
        padding: "10px",
        ...sx
      }}
    >
      <Box sx={{width: { xs: "90%", sm: "90%", md: "60%", lg: "60%", xl: "60%" }, display: "flex", flexDirection: "row", alignItems: "center", gap: "15px"}}>
      {/* <Footbar> */}
        {/* <LinkBox to="/" sx={{ display: "flex", alignItems: "center", gap: "0.7em" }}> */}
          <IconButton
            disableRipple
            sx={{
              bgcolor: "white",
              ":hover": { bgcolor: "primary.dark" }
            }}
          >
            <img width={28} alt="logo" src={logo} />
          </IconButton>
          <Navlink>genomics.ai</Navlink>
        {/* </LinkBox> */}
        {/* <LinkBox to="/about"><Navlink>About us</Navlink></LinkBox>
        <LinkBox to="/docs"><Navlink>Docs </Navlink></LinkBox>
        <LinkBox to="/explore"><Navlink>Explore</Navlink></LinkBox> */}
      {/* </Footbar> */}

      {/* <Box sx={{ width: { xs: "90%", sm: "90%", md: "50%", lg: "50%", xl: "50%" }, height: "2px", bgcolor: "black", margin: "auto" }} /> */}

      <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center" }} >
        <CopyrightIcon sx={{ width: "16px", height: "16px" }} />
        <Typography fontSize="16px" fontWeight="bold" >2022</Typography>
        <Typography fontSize="16px" sx={{ margin: '0 4px' }}> - </Typography>
        <LinkBox to="/imprint" style={{ textDecoration: "none" }}>
          <Typography fontSize="16px" color="black">Imprint</Typography>
        </LinkBox>
        { /*
        <Typography fontSize="16px" sx={{ margin: '0 4px' }}> - </Typography>
        <LinkBox to="/terms" style={{ textDecoration: "none" }}>
          <Typography fontSize="16px" color="black">Terms of Usage</Typography>
        </LinkBox>
        */ }
        <Typography fontSize="16px" sx={{ margin: '0 4px' }}> - </Typography>
        <LinkBox to="/privacy" style={{ textDecoration: "none" }}>
          <Typography fontSize="16px" color="black">Privacy Policy</Typography>
        </LinkBox>
      </Box>
      </Box>
    </Box>
  )
}
