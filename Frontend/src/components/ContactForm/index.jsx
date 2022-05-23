import { Box, Typography, Button, Stack, TextField, TextareaAutosize } from "@mui/material";
import { useState } from "react";
import CustomButton from "components/CustomButton";
import Input from 'components/Input/Input'
import ContactUsService from "shared/services/ContactUs.service";

export default function ContactForm(){

  const [data, setData] = useState({email: "", firstName: "", lastName: "", subject: "", message: ""})
  const [error, setError] = useState({email: false, firstName: false, lastName: false, subject: false, message: false})

  const onChange = (type) => 
    (event) => {
      setData({...data, [type]: event.target.value})
    }

  const onFocus = (type) => 
    () => {
      setError({...error, [type]: false})
    }

  const onBlur = (type) => 
    () => {
    if(data[type]==="") setError({...error, [type]: true})
    else setError({...error, [type]: false})
    }

  const onClick = () => {
    let ok=true
    const newError = {...error}
    for(const type in data){
      if(data[type]===""){
        newError[type]=true
        ok=false
      }
    }
    setError(newError)
    
    ContactUsService.postContactForm(data)
  }

  return (
    <Box sx={{
      width: "100%",
      margin: "2em 0 0em 0",
      padding: "2em 0em",
      boxShadow: "0px 4px 6px rgba(0, 0, 0, 0.10), 0px 0px 1px rgba(0, 0, 0, 0.20)",
      borderRadius: "10px"
    }}>
      <Stack sx={{ width: "80%", margin: "auto" }} direction="column" spacing={4}>
        <Stack direction="row" justifyContent="space-between" spacing={2}>
          <Input isRequired label="Email" errorHandler={error.email} helperText={error.email ? "Email cannot be empty!" : ""} onChangeEvent={onChange("email")} onBlurEvent={onBlur("email")} onFocusEvent={onFocus("email")} />
          <Input isRequired label="First Name" errorHandler={error.firstName} helperText={error.firstName ? "First Name cannot be empty!" : ""} onChangeEvent={onChange("firstName")} onBlurEvent={onBlur("firstName")} onFocusEvent={onFocus("firstName")} />
          <Input isRequired label="Last Name" errorHandler={error.lastName} helperText={error.lastName ? "Last Name cannot be empty!" : ""} onChangeEvent={onChange("lastName")} onBlurEvent={onBlur("lastName")} onFocusEvent={onFocus("lastName")} />
        </Stack>
        <Input isRequired label="Subject" maxLength={100} errorHandler={error.subject} helperText={error.subject ? "Subject cannot be empty!" : ""} onChangeEvent={onChange("subject")} onBlurEvent={onBlur("subject")} onFocusEvent={onFocus("subject")} />
        <Input isRequired label="Message" multiline maxLength={1000} errorHandler={error.message} helperText={error.message ? "Message cannot be empty!" : ""} onChangeEvent={onChange("message")} onBlurEvent={onBlur("message")} onFocusEvent={onFocus("message")} />
        <Box sx={{display: "flex", flexDirection: "row", justifyContent: "center"}}>
          <CustomButton onClick={onClick}>Send</CustomButton>
        </Box>
      </Stack>
    </Box>
  )
}