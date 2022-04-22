import { Box, Typography, Button } from "@mui/material";
import Navbar from "components/Navbar";
import WindowiOS from "components/WindowiOS";
import { useCallback, useState, useEffect, useRef } from "react";
import { colors } from "shared/theme/colors";
import graphic1 from 'assets/landing-illustrations/science.png';
import graphic2 from 'assets/landing-illustrations/upload.png';
import graphic3 from 'assets/landing-illustrations/processing.png';
import graphic4 from 'assets/landing-illustrations/results.png';

const Home = ({ setUser }) => {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false)
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false)

  const onLoginClicked = useCallback(() => {
    setRegistrationFormVisible(false)
    setLoginFormVisible(true)
  }, [setLoginFormVisible])

  const onSignUpClicked = useCallback(() => {
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible])
  
  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  return (
    <Box style={{ overflow: "hidden" }} sx={{ position: "relative" }}>
      {/* STARTING PAGE */}
      <Box sx={{width: "100vw", bgcolor: colors.primary[800], position: "relative", paddingBottom: "4em"}}>
        {/* NAVBAR HERE */}
        <Navbar />
        {/* IOS WINDOW */}
        <WindowiOS />
      </Box>
      {/* the Eclipse */}
      <Box
        sx={{
          position: "relative",
          width: "140vw",
          height: "5vw",
          left: "-20vw", 
          top: "-2vw",
          backgroundColor: "white",
          borderRadius: "50%",
          zIndex: "0"
        }}
      />
      {/* CONTAINER */}
      <Box
        sx={{
          position: "relative",
          width: "100vw",
          paddingLeft: { xs: "10px", md: "15%" },
          paddingRight: { xs: "10px", md: "15%" },
          top: "-2.5vw"
        }}
      >
        
        {/* WHAT WE DO */}
        <Box sx={{ width: "100%" }}>
          <Typography width="100%" sx={{ textAlign: "center" }} fontSize="2em" fontWeight="bold">What we do</Typography>
          <Box sx={{ display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center", width: "80%", margin: "auto" }}>
            <img style={{ width: "25em" }} src={graphic1} alt="Science" />
            <Box sx={{ width: "40%", p: "1em"  }}>
              <Typography fontSize="1.2em" fontWeight="bold">genomics.ai</Typography>
              <Typography margin="2em 0 2em 0">We help you visualize all of your single-cell sequencing data in a fast and easy way with the help of neural networks.</Typography>
              {/* WILL BE REPLACED BY FUTURE BUTTONS */}
              <Button variant="contained">Explore</Button>
            </Box>
          </Box>
        </Box>
        {/* HOW IT WORKS */}
        <Box sx={{ 
          width: "100%", 
          marginTop: "5em", 
          backgroundColor: colors.primary[800],
          p: "1em",
          color: "white",
          borderRadius: "20px",
          boxShadow: "0px 0px 10px rgba(0, 0, 0, 0.1)"
          }}>
          <Typography width="100%" sx={{ textAlign: "center" }} fontSize="2em" fontWeight="bold">How it works</Typography>
          <Box sx={{ gap: "1em", p: "1em", display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center", width: "70%", margin: "auto" }}>
            <Box sx={{ width: "50%" }}>
              <Typography fontSize="1.2em" fontWeight="bold">Upload</Typography>
              <Typography color={colors.neutral[500]}>Upload your files containing the single-cell sequencing data.</Typography>
            </Box>
            <Box sx={{ width: "50%", backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic2} alt="Upload" />
            </Box>
          </Box>  
          <Box sx={{ gap: "1em", p: "1em", display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center", width: "70%", margin: "auto" }}>
            <Box sx={{ width: "50%" }}>
              <Typography fontSize="1.2em" fontWeight="bold">Processing</Typography>
              <Typography color={colors.neutral[500]}>Your input data is processed by our machine learning model which performs the cell-type classification according to your specification.</Typography>
            </Box>
            <Box sx={{ width: "50%", backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic3} alt="Processing" />
            </Box>
          </Box>
          <Box sx={{ gap: "1em", p: "1em", display: "flex", flexDirection: "row", justifyContent: "space-between", alignItems: "center", width: "70%", margin: "auto" }}>
            <Box sx={{ width: "50%" }}>
              <Typography fontSize="1.2em" fontWeight="bold">Check Results</Typography>
              <Typography color={colors.neutral[500]}>Your input data is processed by our machine learning model which performs the cell-type classification according to your specification.After the algorithm has processed the data you can specify the project you want your cell-type data to be associated with and view the results.</Typography>
            </Box>
            <Box sx={{ width: "50%", backgroundColor: "white", borderRadius: "20px" }}>
              <img style={{ width: "100%" }} src={graphic4} alt="Check Results" />
            </Box>
          </Box>
        </Box>
      </Box>

      {/* FOOTER */}
      
    </Box>
  )
}

export default Home