import * as React from 'react';
import { styled } from '@mui/material/styles';
import Button from '@mui/material/Button';
import UploadFileIcon from '@mui/icons-material/UploadFile';

const Input = styled('input')({
  display: 'none',
});

function UploadButton({ onChange, disabled }) {
  let button;

  // TODO: Looks sus
  if (!disabled) {
    button = (
      <Button component="span" endIcon={<UploadFileIcon />}>
        Choose a file
      </Button>
    );
  } else {
    button = (
      <Button component="span" disabled={disabled} endIcon={<UploadFileIcon />}>
        Choose a file
      </Button>
    );
  }

  return (
    <label htmlFor="contained-button-file">
      <Input
        accept=".h5ad"
        id="contained-button-file"
        multiple
        type="file"
        onChange={onChange}
        disabled={disabled}
      />
      {button}
    </label>
  );
}

export default UploadButton;
