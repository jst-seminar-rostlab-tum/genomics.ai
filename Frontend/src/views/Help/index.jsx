import React, { useCallback, useState } from 'react';
import { Box, Typography, Button, Stack, TextField, TextareaAutosize } from "@mui/material";
import Input from 'components/Input/Input'
import HeaderView from 'components/general/HeaderView';
import CustomButton from "components/CustomButton";
import ContactForm from 'components/ContactForm';

export default function Help({ sidebarShown }) {
  const [contactDetails, setContactDetails] = useState({
    message: '',
  });

  const handleTextChange = useCallback((e) => {
    setContactDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [setContactDetails]);

  const formSubmit = async (e) => {
    e.preventDefault();

    const data = {
      message: contactDetails.message,
    };

    try {
      // send POST to the server
      // await axios.post("BACKEND_URL", data);
      setContactDetails({
        message: '...sending',
      });
    } catch (error) {
      console.log(error);
    }
  };

  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title="Contact Us"
    >

      <Box sx={{margin: "1em auto", display: "flex", flexDirection: "column", alignItems: "center", textAlign: "center"}}>
        <Typography fontWeight={400} fontSize="1.2em">Please describe to us the problem you are experiencing in the contact form below.</Typography>
      </Box>

      <Box sx={{margin: "1em auto", display: "flex", flexDirection: "column", alignItems: "center", textAlign: "center"}}>
        <Typography fontSize="1em">We will get back to you as soon as possible.</Typography>
      </Box>

      <Box sx={{ position: "relative", width: "60vw", margin: "1em auto" }} >
        {/* <Box sx={{
          width: "100%",
          margin: "2em 0 0em 0",
          padding: "2em 0em",
          boxShadow: "0px 4px 6px rgba(0, 0, 0, 0.10), 0px 0px 1px rgba(0, 0, 0, 0.20)",
          borderRadius: "10px"
        }}>
          <Stack sx={{ width: "80%", margin: "auto" }} direction="column" spacing={4}>
            <Stack direction="row" justifyContent="space-between" spacing={2}>
              <Input isRequired label="Email" />
              <Input isRequired label="First Name" />
              <Input isRequired label="Last Name" />
            </Stack>
            <Input label="Subject" maxLength={100} />
            <Input isRequired label="Message" multiline maxLength={1000} />
            <Box sx={{display: "flex", flexDirection: "row", justifyContent: "center"}}>
              <CustomButton>Send</CustomButton>
            </Box>
          </Stack>
        </Box> */}
        <ContactForm />
      </Box>
    </HeaderView>
  );
}
