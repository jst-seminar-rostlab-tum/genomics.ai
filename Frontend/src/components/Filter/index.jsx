import React from 'react';
import FormGroup from '@mui/material/FormGroup';
import Box from '@mui/material/Box';
import MenuItem from '@mui/material/MenuItem';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import FormControlLabel from '@mui/material/FormControlLabel';
import Checkbox from '@mui/material/Checkbox';
import styles from './filter.module.css';
import { Divider, Stack, TextField } from '@mui/material';

export default function Filter() {
  const [reference, setReference] = React.useState('');

  const handleReferenceChange = (event) => {
    setReference(event.target.value);
  };

  const [category, setCategory] = React.useState('');

  const handleCategoryChange = (event) => {
    setCategory(event.target.value);
  };

  return (
    <Box sx={{ width: 420 }}>
      <FormGroup>
        <Stack
          className="flexContainer"
          direction="row"
          sx={{ margin: 4 }}

        >
          <FormControlLabel control={<Checkbox />} label={<Typography sx={{ fontWeight: 'bold' }}>Referenced By Atlas</Typography>} />
          <TextField
            autoWidth
            value={reference}
            onChange={handleReferenceChange}
            sx={{ minWidth: '160px' }}
            select
            variant="standard"
          >
            <MenuItem value="Human - PBMC">Human - PBMC</MenuItem>
            <MenuItem value="Human - PBMC2">Human - PBMC2</MenuItem>
            <MenuItem value="Human - PBMC3">Human - PBMC3</MenuItem>
          </TextField>
        </Stack>
        <Divider />
        <Stack
          className="flexContainer"
          direction="row"
          sx={{ margin: 4 }}
        >
          <FormControlLabel control={<Checkbox />} label={<Typography sx={{ fontWeight: 'bold' }}>Category</Typography>} />
          <TextField
            autoWidth
            value={category}
            onChange={handleCategoryChange}
            sx={{ minWidth: '160px' }}
            select
            variant="standard"
          >
            <MenuItem value="lung">lung</MenuItem>
            <MenuItem value="lung2">lung2</MenuItem>
            <MenuItem value="lung3">lung3</MenuItem>
          </TextField>
          <Button variant="contained">Add</Button>
        </Stack>
      </FormGroup>

    </Box>
  );
}
