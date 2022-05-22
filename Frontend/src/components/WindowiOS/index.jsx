/* eslint-disable */

import { Circle } from "@mui/icons-material"
import { Box, Typography } from "@mui/material"
import { colors } from "shared/theme/colors"
import tum from 'assets/landing-illustrations/tum-logo.png';
import rostlab from 'assets/landing-illustrations/rostlab.png';
import helmholtz from 'assets/landing-illustrations/helmholtz.png';
import CustomButton from "components/CustomButton";
import { useHistory } from "react-router-dom";

const WindowiOS = ({ onSignUpClicked }) => {
  const history = useHistory()

  return (
    <Box
      sx={{
        position: "relative",
        width: { xs: "90%", sm: "90%", md: "61.8%", lg: "61.8%", xl: "61.8%" },
        margin: "auto",
        backgroundColor: "rgba(20, 20, 20, 0.6)",
        color: "white",
        p: "0.7em",
        borderRadius: "20px",
        boxShadow: "0px 0px 10px rgba(255,255,255, 0.15)",
        marginTop: "3%",
        marginBottom: "2em"
      }}
    >
      <Box
        sx={{
          position: "absolute",
          display: "flex",
          flexDirection: "row",
          gap: "10px",
          left: "3%",
        }}
      >
        <Circle sx={{ width: "15px", color: colors.error.main }} />
        <Circle sx={{ width: "15px", color: colors.warning.main }} />
        <Circle sx={{ width: "15px", color: colors.success.main }} />
      </Box>
      <Typography textAlign="center">genomics.ai</Typography>
      <Box sx={{ p: { xs: "1em 1em", sm: "1em 4em", md: "1em 4em", lg: "1em 4em", xl: "1em 4em" }, marginTop: "4em" }}>
        <Typography fontSize={{ xs: "1.7rem", sm: "2.3rem", md: "2.1rem", lg: "3.1rem", xl: "3.1rem" }} fontWeight="bold">
          AI-Driven Cell Type Annotation
        </Typography>
        <Typography fontSize={{ xs: "1.4rem", sm: "1.6rem", md: "1.4rem", lg: "2rem", xl: "2rem" }} fontWeight="semibold">
          No Code. Just Results.
        </Typography>
        <Typography fontSize="1.2rem" fontWeight="light" color="#8193B2">
          Genomics.ai helps you visualize all of your single-cell sequencing data in a fast and easy way using neural networks.
        </Typography>
        <Box sx={{ margin: "2em 0 0 0px", display: "flex", flexDirection: "row", gap: "10px", alignItems: "center" }}>
          <CustomButton onClick={onSignUpClicked} sx={{ width: "150px" }}>Sign up</CustomButton>
        </Box>
        <Box sx={{ marginBlock: "10px", display: "flex", flexDirection: "row", gap: "10px", alignItems: "center" }}>
          <Typography sx={{ color: colors.neutral[200] }}>Test out our GeneMapper. No need to sign up.</Typography>
          <CustomButton type="secondary" onClick={() => history.push("/explore")} sx={{ color: colors.secondary1[700] }}>Explore</CustomButton>
        </Box>
        
        {/* DIVIDER */}
        <Box sx={{ marginTop: "3em", display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center" }}>
          <Box sx={{ height: "1px", width: { xs: "25%", sm: "37%", md: "35%", lg: "40%", xl: "40%" }, backgroundColor: "rgba(255,255,255, 0.5)" }} />
          <Typography fontSize={{ xs: "0.6em", sm: "1em", md: "1em", lg: "1em", xl: "1em" }} >Our partners</Typography>
          <Box sx={{ height: "1px", width: { xs: "25%", sm: "37%", md: "35%", lg: "40%", xl: "40%" }, backgroundColor: "rgba(255,255,255, 0.5)" }} />
        </Box>
      </Box>
      {/* BUTTON PLACEHOLDER UNTIL WE HAVE A BUTTON DESIGN */}
      <Box
        sx={{
          display: "flex",
          flexDirection: { xs: "column", sm: "row", md: "row", lg: "row", xl: "row" },
          justifyContent: "space-between",
          alignItems: "center",
          width: "80%",
          margin: "auto",
          p: "1em 1em 2em 1em"
        }}
      >
        <Box sx={{ margin: { xs: "2em 2em", sm: "0em", md: "0em", lg: "0em", xl: "0em" } }} >
          <img src={helmholtz} alt="HELMHOLTZ" style={{ width: "8em", height: "auto" }} />
        </Box>
        <Box sx={{ margin: { xs: "2em 2em", sm: "0em", md: "0em", lg: "0em", xl: "0em" } }} >
          <img src={rostlab} alt="ROSTLAB" style={{ width: "10em", height: "auto" }} />
        </Box>
        <Box sx={{ margin: { xs: "2em 2em", sm: "0em", md: "0em", lg: "0em", xl: "0em" } }} >
          <img src={tum} alt="TUM" style={{ width: "5em", height: "auto" }} />
        </Box>
      </Box>
    </Box>
  )
}

export default WindowiOS