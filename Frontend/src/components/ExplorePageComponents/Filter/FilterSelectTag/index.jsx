import React, { useState } from 'react';
import {
  Stack, Button, Box, TextField, MenuItem, Typography,
} from '@mui/material';
import Tag from 'components/Tag';

// General filter that is needed in all categories
export default function FilterSelectTag(props) {
  const [category, setCategory] = useState('');
  const [selectedCategories, setSelectedCategories] = useState([]);
  const { categories, onChange, label } = props;

  const handleCategoryChange = (event) => {
    setCategory(event.target.value);
  };

  const handleCategoriesAppend = (newValue) => {
    setSelectedCategories((array) => [...array, newValue]);
    onChange(selectedCategories);
  };

  const addCategory = () => {
    if (!selectedCategories.includes(category) && category !== '') {
      handleCategoriesAppend(category);
    }
  };

  const removeTag = (selectedCategory) => {
    setSelectedCategories(selectedCategories.filter((item) => item !== selectedCategory));
    onChange(selectedCategories);
  };
  return (
    <Box>
      <Stack
        direction="row"
        spacing={2}
        alignItems="center"
        justifyContent="flex-end"
        sx={{ marginBottom: 2 }}
      >
        <Typography sx={{ fontWeight: 'bold' }}>{label}</Typography>
        <TextField
          value={category}
          onChange={handleCategoryChange}
          select
          variant="standard"
        >
          {categories.map((categoryItem) => (
            <MenuItem key={categoryItem} value={categoryItem}>{categoryItem}</MenuItem>
          ))}
        </TextField>
        <Button variant="contained" onClick={addCategory}>Add</Button>
      </Stack>
      <Stack
        direction="column"
        spacing={2}
      >
        {selectedCategories.map((sCategory, index) => (
          <Tag
            key={index}
            content={sCategory}
            variant="primary-default"
            isDeletable
            handleDeletion={() => removeTag(sCategory)}
          />
        ))}
      </Stack>
    </Box>
  );
}
