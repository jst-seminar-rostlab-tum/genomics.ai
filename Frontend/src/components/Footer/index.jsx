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
        height: "40px",
        display: "flex",
        flexDirection: "row",
        alignItems: "center",
        justifyContent: "center",
        padding: "10px",
        ...sx
      }}
    >
      <Box sx={{display: "flex", flexDirection: "row", alignItems: "center", gap: "15px", alignSelf: "center"}}>
      <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center" }} >
        <Typography fontSize="12px" fontWeight="bold" >© Copyright 2022</Typography>
        <Typography fontSize="12px" sx={{ margin: '0 4px' }}> - </Typography>
        <LinkBox to="/legalnotice" style={{ textDecoration: "none" }}>
          <Typography fontSize="12px" color="black">Legal Notice</Typography>
        </LinkBox>
        { /*
        <Typography fontSize="16px" sx={{ margin: '0 4px' }}> - </Typography>
        <LinkBox to="/terms" style={{ textDecoration: "none" }}>
          <Typography fontSize="16px" color="black">Terms of Usage</Typography>
        </LinkBox>
        */ }
        <Typography fontSize="12px" sx={{ margin: '0 4px' }}> - </Typography>
        <LinkBox to="/privacy" style={{ textDecoration: "none" }}>
          <Typography fontSize="12px" color="black">Privacy Policy</Typography>
        </LinkBox>
      </Box>
      </Box>
    </Box>
  )
}
