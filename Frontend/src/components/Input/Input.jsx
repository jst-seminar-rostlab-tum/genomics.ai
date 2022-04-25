import * as React from 'react';
import Box from '@mui/material/Box';
import TextField from '@mui/material/TextField';

/* Props: Most of MUI TextField props..
Example usage
Import Input to any component, see below examples:

const [name, setName] = useState('');

<Input isRequired onChangeEvent={setName} disabledHandler={false} />
<Input isRequired onChangeEvent={setName} disabledHandler={false} />
<Input helperText="Incorrect Output" errorHandler isRequired
 onChangeEvent={setName} disabledHandler={false} />

 //For Multiline
 <Input isRequired maxLength={120} multiline disabledHandler={false} />
 */

const Input = ({
  helperText, placeholder, errorHandler, disabledHandler,
  isRequired = true, defaultValue, multiline = false, onChangeEvent = null,
  label = 'Required',
  maxLength = 40,
  type = 'text',
}) => (
  <Box
    component="form"
    sx={{
      '& .MuiTextField-root': { m: 1, width: '100%' },
      '& .MuiFormLabel-asterisk': { color: '#FF5F58' },
    }}
    noValidate
    autoComplete="off"
  >
    <TextField
      error={errorHandler}
      required={isRequired}
      disabled={disabledHandler}
      placeholder={placeholder}
      defaultValue={defaultValue}
      inputProps={{ maxLength }}
      multiline={multiline}
      type={type}
      id="input"
      label={label}
      helperText={helperText}
      variant="standard"
      onChange={(e) => (onChangeEvent !== null ? onChangeEvent(e.target.value) : null)}
    />
  </Box>
);

export default Input;
